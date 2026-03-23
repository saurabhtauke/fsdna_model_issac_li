#!/usr/bin/env python3
"""
Predict the extension of a flow-stretched dsDNA tether in a microfluidic channel.

This script implements the model described in:
"Massively Parallel Bead-Free Force Spectroscopy with Fluorescence"
(bioRxiv 2025.09.14.675960), specifically the "Modelling Flow-Stretched DNA"
section.

How to use
----------
1. Edit the file `flow_stretched_dna_params.py` in the same directory.
2. Run:
       python flow_stretched_dna_model.py

No command-line inputs are required.

Model summary
-------------
1. Convert flow rate Q to wall shear stress using
       sigma = 6 * eta * Q / (h^2 * w)
   for a rectangular channel.
2. Discretize the DNA contour into ~150 bp segments.
3. Assign each segment a height profile
       h_i = alpha * sigma^beta * d_i^gamma
   where d_i is the segment-center contour position measured from the free end,
   as a fraction of total contour length.
4. Compute segment drag
       f_i = sigma * n_i * h_i
   using sigma in Pa and n_i, h_i in microns, which returns f_i in pN because
   1 Pa * um^2 = 1 pN.
5. Compute the tension profile from the free end toward the anchor.
6. Convert the average tension on each segment into an axial extension using an
   extensible worm-like chain (extWLC / XWLC) model.
7. Sum segment extensions to obtain total tether extension.

XWLC implementation
-------------------
The force-extension inversion is implemented using the same extWLC functional
form used in `wlc/models.py` from:
https://github.com/franciscopalmeromoya/wlc-fitting

    F = (kBT / Lp) * [1 / (4 (1 - d/Lc + F/S)^2) - 1/4 + d/Lc - F/S]

Here the script solves that relation numerically for d at a given segment force
F, so the script stays self-contained and does not require lmfit.

Important implementation note
-----------------------------
The paper excerpt states that segment extension is obtained numerically from the
inverse XWLC relation, but it does not specify which within-segment tension is
used when tension varies across a segment. Here I use the average of the
segment's upstream and downstream tensions, which is the physically natural
choice for a finite segment.
"""

from __future__ import annotations

import csv
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import List, Mapping, Sequence


# Physical constants
K_B = 1.380649e-23  # J / K
NM_PER_BP = 0.34    # contour rise of B-DNA in nm / bp
UM_PER_BP = NM_PER_BP * 1e-3


@dataclass
class ModelParams:
    # Paper-fit parameters
    alpha: float = 0.160
    beta: float = -0.229
    gamma: float = -0.157
    nominal_segment_bp: int = 150
    temperature_c: float = 20.0

    # Looked-up / user-overridable DNA mechanics parameters
    persistence_length_nm: float = 45.0
    stretch_modulus_pN: float = 1100.0

    # Fluid property
    viscosity_pa_s: float = 1.002e-3

    @property
    def temperature_K(self) -> float:
        return self.temperature_c + 273.15

    @property
    def kbt_pN_nm(self) -> float:
        return K_B * self.temperature_K * 1e21  # J -> pN*nm


@dataclass
class SegmentState:
    index_from_free_end: int
    bp: int
    contour_um: float
    d_fraction: float
    height_um: float
    drag_pN: float
    tension_free_side_pN: float
    tension_anchor_side_pN: float
    tension_avg_pN: float
    extension_um: float


@dataclass
class Prediction:
    flow_rate_ul_min: float
    shear_stress_pa: float
    anchor_force_pN: float
    extension_um: float
    contour_length_um: float
    fractional_extension: float
    n_segments: int
    segments: List[SegmentState]


@dataclass
class RunConfig:
    flow_rates_ul_min: List[float]
    channel_height_um: float
    channel_width_um: float
    dna_length_bp: int
    params: ModelParams
    csv_out: str | None = None
    dump_segments_for_first_flow: bool = False


def ul_min_to_m3_s(q_ul_min: float) -> float:
    return q_ul_min * 1e-9 / 60.0


def um_to_m(x_um: float) -> float:
    return x_um * 1e-6


def contour_length_um_from_bp(bp: int) -> float:
    return bp * UM_PER_BP


def discretize_dna_bp(total_bp: int, nominal_segment_bp: int) -> List[int]:
    if total_bp <= 0:
        raise ValueError("DNA length must be positive.")
    if nominal_segment_bp <= 0:
        raise ValueError("nominal_segment_bp must be positive.")

    if total_bp <= nominal_segment_bp:
        return [total_bp]

    n_full = total_bp // nominal_segment_bp
    remainder = total_bp % nominal_segment_bp

    if remainder == 0:
        return [nominal_segment_bp] * n_full

    if n_full == 1:
        return [total_bp]

    # The paper says any remainder is assigned to the segment closest to the anchor.
    return [nominal_segment_bp] * (n_full - 1) + [nominal_segment_bp + remainder]



def wall_shear_stress_pa(
    flow_rate_ul_min: float,
    channel_height_um: float,
    channel_width_um: float,
    viscosity_pa_s: float,
) -> float:
    q_m3_s = ul_min_to_m3_s(flow_rate_ul_min)
    h_m = um_to_m(channel_height_um)
    w_m = um_to_m(channel_width_um)
    if h_m <= 0 or w_m <= 0:
        raise ValueError("Channel dimensions must be positive.")
    return 6.0 * viscosity_pa_s * q_m3_s / (h_m * h_m * w_m)



def extwlc_rhs_pN(
    force_pN: float,
    extension_um: float,
    contour_um: float,
    persistence_length_nm: float,
    stretch_modulus_pN: float,
    kbt_pN_nm: float,
) -> float:
    """
    extWLC force law mirroring the relation in wlc-fitting/wlc/models.py.

    The repository implementation uses distance in um and contour length in nm;
    this helper keeps the same physics but accepts the segment contour length in
    um for convenience.
    """
    if contour_um <= 0.0:
        raise ValueError("Contour length must be positive.")
    if persistence_length_nm <= 0.0 or stretch_modulus_pN <= 0.0:
        raise ValueError("Persistence length and stretch modulus must be positive.")

    d_nm = extension_um * 1000.0
    lc_nm = contour_um * 1000.0
    l = d_nm / lc_nm - force_pN / stretch_modulus_pN

    return (kbt_pN_nm / persistence_length_nm) * (
        0.25 / (1.0 - l) ** 2 - 0.25 + l
    )



def extwlc_extension_fraction(
    force_pN: float,
    contour_um: float,
    persistence_length_nm: float,
    stretch_modulus_pN: float,
    kbt_pN_nm: float,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> float:
    """
    Numerically invert the extWLC relation to obtain x/Lc for a given force.

    The bracketing interval follows directly from the extWLC normalized variable
        l = d/Lc - F/S
    with the physical singular limit l -> 1^-.
    """
    if force_pN < 0.0:
        raise ValueError("Force must be non-negative.")
    if contour_um <= 0.0:
        raise ValueError("Contour length must be positive.")
    if persistence_length_nm <= 0.0 or stretch_modulus_pN <= 0.0:
        raise ValueError("Persistence length and stretch modulus must be positive.")

    if force_pN == 0.0:
        return 0.0

    def residual(extension_um: float) -> float:
        return force_pN - extwlc_rhs_pN(
            force_pN=force_pN,
            extension_um=extension_um,
            contour_um=contour_um,
            persistence_length_nm=persistence_length_nm,
            stretch_modulus_pN=stretch_modulus_pN,
            kbt_pN_nm=kbt_pN_nm,
        )

    lo = 0.0
    singular_limit_um = contour_um * (1.0 + force_pN / stretch_modulus_pN)
    hi = singular_limit_um * (1.0 - 1e-12)

    r_lo = residual(lo)
    r_hi = residual(hi)
    if r_lo < 0.0 or r_hi > 0.0:
        raise RuntimeError("Failed to bracket extWLC inversion root.")

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        r_mid = residual(mid)
        if abs(r_mid) < tol or (hi - lo) < tol:
            return mid / contour_um
        if r_mid > 0.0:
            lo = mid
        else:
            hi = mid

    return 0.5 * (lo + hi) / contour_um



def predict_extension(
    flow_rate_ul_min: float,
    channel_height_um: float,
    channel_width_um: float,
    dna_length_bp: int,
    params: ModelParams,
) -> Prediction:
    sigma_pa = wall_shear_stress_pa(
        flow_rate_ul_min=flow_rate_ul_min,
        channel_height_um=channel_height_um,
        channel_width_um=channel_width_um,
        viscosity_pa_s=params.viscosity_pa_s,
    )

    contour_total_um = contour_length_um_from_bp(dna_length_bp)
    segment_bps = discretize_dna_bp(dna_length_bp, params.nominal_segment_bp)

    if sigma_pa <= 0.0:
        return Prediction(
            flow_rate_ul_min=flow_rate_ul_min,
            shear_stress_pa=0.0,
            anchor_force_pN=0.0,
            extension_um=0.0,
            contour_length_um=contour_total_um,
            fractional_extension=0.0,
            n_segments=len(segment_bps),
            segments=[],
        )

    segments: List[SegmentState] = []
    cumulative_contour_um = 0.0
    cumulative_drag_pN = 0.0
    total_extension_um = 0.0

    for idx, bp in enumerate(segment_bps, start=1):
        contour_um = contour_length_um_from_bp(bp)
        center_um = cumulative_contour_um + 0.5 * contour_um
        d_fraction = center_um / contour_total_um

        height_um = params.alpha * (sigma_pa ** params.beta) * (d_fraction ** params.gamma)
        drag_pN = sigma_pa * contour_um * height_um

        tension_free_side_pN = cumulative_drag_pN
        tension_anchor_side_pN = cumulative_drag_pN + drag_pN
        tension_avg_pN = 0.5 * (tension_free_side_pN + tension_anchor_side_pN)

        extension_fraction = extwlc_extension_fraction(
            force_pN=tension_avg_pN,
            contour_um=contour_um,
            persistence_length_nm=params.persistence_length_nm,
            stretch_modulus_pN=params.stretch_modulus_pN,
            kbt_pN_nm=params.kbt_pN_nm,
        )
        extension_um = contour_um * extension_fraction

        segments.append(
            SegmentState(
                index_from_free_end=idx,
                bp=bp,
                contour_um=contour_um,
                d_fraction=d_fraction,
                height_um=height_um,
                drag_pN=drag_pN,
                tension_free_side_pN=tension_free_side_pN,
                tension_anchor_side_pN=tension_anchor_side_pN,
                tension_avg_pN=tension_avg_pN,
                extension_um=extension_um,
            )
        )

        cumulative_contour_um += contour_um
        cumulative_drag_pN += drag_pN
        total_extension_um += extension_um

    anchor_force_pN = cumulative_drag_pN
    fractional_extension = total_extension_um / contour_total_um if contour_total_um > 0 else 0.0

    return Prediction(
        flow_rate_ul_min=flow_rate_ul_min,
        shear_stress_pa=sigma_pa,
        anchor_force_pN=anchor_force_pN,
        extension_um=total_extension_um,
        contour_length_um=contour_total_um,
        fractional_extension=fractional_extension,
        n_segments=len(segment_bps),
        segments=segments,
    )



def format_predictions_table(predictions: Sequence[Prediction]) -> str:
    headers = [
        "flow_uL_min",
        "sigma_Pa",
        "anchor_force_pN",
        "extension_um",
        "contour_um",
        "frac_ext",
        "n_segments",
    ]

    rows = [
        [
            f"{p.flow_rate_ul_min:.6g}",
            f"{p.shear_stress_pa:.6g}",
            f"{p.anchor_force_pN:.6g}",
            f"{p.extension_um:.6g}",
            f"{p.contour_length_um:.6g}",
            f"{p.fractional_extension:.6g}",
            str(p.n_segments),
        ]
        for p in predictions
    ]

    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]

    def fmt_row(row: Sequence[str]) -> str:
        return "  ".join(cell.rjust(widths[i]) for i, cell in enumerate(row))

    lines = [fmt_row(headers), fmt_row(["-" * w for w in widths])]
    lines.extend(fmt_row(row) for row in rows)
    return "\n".join(lines)



def write_csv(predictions: Sequence[Prediction], csv_path: str) -> None:
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "flow_rate_uL_min",
                "shear_stress_Pa",
                "anchor_force_pN",
                "extension_um",
                "contour_length_um",
                "fractional_extension",
                "n_segments",
            ]
        )
        for p in predictions:
            writer.writerow(
                [
                    p.flow_rate_ul_min,
                    p.shear_stress_pa,
                    p.anchor_force_pN,
                    p.extension_um,
                    p.contour_length_um,
                    p.fractional_extension,
                    p.n_segments,
                ]
            )



def print_segment_diagnostics(prediction: Prediction) -> None:
    headers = [
        "seg",
        "bp",
        "d_frac",
        "height_um",
        "drag_pN",
        "F_free_pN",
        "F_anchor_pN",
        "F_avg_pN",
        "ext_um",
    ]
    rows = [
        [
            str(s.index_from_free_end),
            str(s.bp),
            f"{s.d_fraction:.6g}",
            f"{s.height_um:.6g}",
            f"{s.drag_pN:.6g}",
            f"{s.tension_free_side_pN:.6g}",
            f"{s.tension_anchor_side_pN:.6g}",
            f"{s.tension_avg_pN:.6g}",
            f"{s.extension_um:.6g}",
        ]
        for s in prediction.segments
    ]
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]

    def fmt_row(row: Sequence[str]) -> str:
        return "  ".join(cell.rjust(widths[i]) for i, cell in enumerate(row))

    print("\nSegment diagnostics (first flow rate only):")
    print(fmt_row(headers))
    print(fmt_row(["-" * w for w in widths]))
    for row in rows:
        print(fmt_row(row))



def _load_module_from_path(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module



def _get_required(config: Mapping[str, object], key: str):
    if key not in config:
        raise KeyError(f"Missing required parameter '{key}' in parameters file.")
    return config[key]



def load_run_config(params_file: Path) -> RunConfig:
    if not params_file.exists():
        raise FileNotFoundError(
            f"Parameters file not found: {params_file}\n"
            "Create or edit flow_stretched_dna_params.py and try again."
        )

    module = _load_module_from_path("flow_stretched_dna_params", params_file)

    if not hasattr(module, "CONFIG"):
        raise AttributeError(
            f"Parameters file {params_file} must define a CONFIG dictionary."
        )

    config = getattr(module, "CONFIG")
    if not isinstance(config, dict):
        raise TypeError("CONFIG must be a Python dictionary.")

    flow_rates_ul_min = list(_get_required(config, "flow_rates_ul_min"))
    channel_height_um = float(_get_required(config, "channel_height_um"))
    channel_width_um = float(_get_required(config, "channel_width_um"))
    dna_length_bp = int(_get_required(config, "dna_length_bp"))

    params = ModelParams(
        alpha=float(config.get("alpha", 0.160)),
        beta=float(config.get("beta", -0.229)),
        gamma=float(config.get("gamma", -0.157)),
        nominal_segment_bp=int(config.get("nominal_segment_bp", 150)),
        temperature_c=float(config.get("temperature_c", 20.0)),
        persistence_length_nm=float(config.get("persistence_length_nm", 45.0)),
        stretch_modulus_pN=float(config.get("stretch_modulus_pN", 1100.0)),
        viscosity_pa_s=float(config.get("viscosity_pa_s", 1.002e-3)),
    )

    csv_out = config.get("csv_out")
    if csv_out is not None:
        csv_out = str(csv_out)

    return RunConfig(
        flow_rates_ul_min=[float(x) for x in flow_rates_ul_min],
        channel_height_um=channel_height_um,
        channel_width_um=channel_width_um,
        dna_length_bp=dna_length_bp,
        params=params,
        csv_out=csv_out,
        dump_segments_for_first_flow=bool(config.get("dump_segments_for_first_flow", False)),
    )



def main() -> None:
    script_dir = Path(__file__).resolve().parent
    params_path = script_dir / "flow_stretched_dna_params.py"
    run_config = load_run_config(params_path)

    predictions = [
        predict_extension(
            flow_rate_ul_min=flow,
            channel_height_um=run_config.channel_height_um,
            channel_width_um=run_config.channel_width_um,
            dna_length_bp=run_config.dna_length_bp,
            params=run_config.params,
        )
        for flow in run_config.flow_rates_ul_min
    ]

    print(f"Loaded parameters from: {params_path}")
    print(format_predictions_table(predictions))

    if run_config.dump_segments_for_first_flow and predictions:
        print_segment_diagnostics(predictions[0])

    if run_config.csv_out:
        csv_path = Path(run_config.csv_out)
        if not csv_path.is_absolute():
            csv_path = script_dir / csv_path
        write_csv(predictions, str(csv_path))
        print(f"\nWrote CSV summary to: {csv_path}")


if __name__ == "__main__":
    main()

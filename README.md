# fsdna-model-issac-li

Flow-stretched dsDNA extension model for a tethered molecule in a rectangular microfluidic channel.

## What this model does

Given flow rate, channel geometry, and DNA length, the model predicts:
- wall shear stress
- per-segment drag and tension profile
- anchor force (pN)
- total DNA extension (um)
- fractional extension (extension / contour length)

The implementation follows the "Modelling Flow-Stretched DNA" section described in:
"Massively Parallel Bead-Free Force Spectroscopy with Fluorescence" (bioRxiv 2025.09.14.675960).

## Model summary

For each flow condition, the script computes:

1. Wall shear stress for a rectangular channel:

   `sigma = 6 * eta * Q / (h^2 * w)`

2. DNA discretization into nominal 150 bp segments.

3. Segment height profile:

   `h_i = alpha * sigma^beta * d_i^gamma`

   where `d_i` is segment-center contour position from the free end, normalized by contour length.

4. Segment drag force:

   `f_i = sigma * n_i * h_i`

   using `sigma` in Pa and `n_i`, `h_i` in um, giving force directly in pN (`1 Pa * um^2 = 1 pN`).

5. Tension accumulation from free end to anchor.

6. Segment extension from an extensible WLC (extWLC/XWLC) inversion using segment-average tension.

7. Total extension by summing segment extensions.

## Repository files

- `flow_stretched_dna_model.py`: main model and prediction pipeline
- `flow_stretched_dna_params.py`: editable `CONFIG` inputs
- `flow_stretched_dna_explorer.ipynb`: notebook workflow
- `flow_stretched_dna_notes.md`: short notes
- `prediction.csv`: example output

## Inputs (CONFIG)

Required:
- `flow_rates_ul_min` (list of uL/min values)
- `channel_height_um`
- `channel_width_um`
- `dna_length_bp`

Common optional:
- `csv_out` (filename or `None`)
- `dump_segments_for_first_flow` (`True`/`False`)
- `alpha`, `beta`, `gamma`
- `nominal_segment_bp`
- `temperature_c`
- `persistence_length_nm`
- `stretch_modulus_pN`
- `viscosity_pa_s`

## Install and run

Using PDM:

```bash
pdm install
pdm run python flow_stretched_dna_model.py
```

Using an existing virtual environment:

```bash
source .venv/bin/activate
python flow_stretched_dna_model.py
```

## Output

The script prints a per-flow summary table with:
- `flow_uL_min`
- `sigma_Pa`
- `anchor_force_pN`
- `extension_um`
- `contour_um`
- `frac_ext`
- `n_segments`

If `csv_out` is set, it also writes the same summary to CSV.

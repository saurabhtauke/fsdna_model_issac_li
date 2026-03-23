# Flow-stretched DNA model notes

## Files
- `flow_stretched_dna_model.py`: main model script.
- `flow_stretched_dna_params.py`: editable parameter file.

## Usage
1. Edit `flow_stretched_dna_params.py`.
2. Run:
   ```bash
   python flow_stretched_dna_model.py
   ```

## Required inputs in `CONFIG`
- `flow_rates_ul_min`: one or more flow rates in uL/min.
- `channel_height_um`: channel height in um.
- `channel_width_um`: channel width in um.
- `dna_length_bp`: DNA length in base pairs.

## Common optional fields
- `csv_out`: output CSV filename, or `None` to disable CSV export.
- `dump_segments_for_first_flow`: print per-segment diagnostics for the first flow rate.
- `persistence_length_nm`: default 45 nm.
- `stretch_modulus_pN`: default 1100 pN.
- `viscosity_pa_s`: default 1.002e-3 Pa s.

## Model structure
The script keeps the paper-derived hydrodynamic model and uses a numerical inversion of the extWLC relation for segment extension.

"""Editable parameters for flow_stretched_dna_model.py.

Edit the numbers below, save the file, then run:
    python flow_stretched_dna_model.py
"""

CONFIG = {
    # Required experimental inputs
    "flow_rates_ul_min": [0.1, 0.5, 1.0, 2.0, 5.0, 10.0],
    "channel_height_um": 120.0,
    "channel_width_um": 1200.0,
    "dna_length_bp": 22900,

    # Optional output controls
    "csv_out": "prediction.csv",   # set to None to disable CSV output
    "dump_segments_for_first_flow": False,

    # Paper-fit parameters
    "alpha": 0.160,
    "beta": -0.229,
    "gamma": -0.157,
    "nominal_segment_bp": 150,
    "temperature_c": 21.5,

    # DNA mechanics parameters
    "persistence_length_nm": 45.0,
    "stretch_modulus_pN": 1100.0,

    # Fluid property
    "viscosity_pa_s": 1.002e-3,
}

# Shared YAML-reader defaults. Export cold metadata with tools/export-node-defaults.jl.
const _NODE_DEFAULTS = (
    cmos=(;binning=1, column_readout_noise_e=0.0, row_readout_noise_e=0.0, photon_noise=true),
    emccd=(;photon_noise=true, dark_current_e_per_pixel_s=0),
    pyramid=(;modulation_propagation_strategy="pupil_tilt",
        pupil_shift_x_pixels=(0, 0, 0, 0), pupil_shift_y_pixels=(0, 0, 0, 0)),
)

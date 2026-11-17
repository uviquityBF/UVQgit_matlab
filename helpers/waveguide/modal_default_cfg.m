function cfg = modal_default_cfg(cfg_override)
% Default settings shared by modal_power_distribution / modal_power_sweep.
% Any field of the optional struct cfg_override replaces the default.
%
% USAGE:
%   cfg = modal_default_cfg();
%   cfg = modal_default_cfg(struct('mode_file', 'C:\...\modes.csv'));

    cfg.spot_source        = 'measured'; % 'zemax' | 'gaussian' | 'measured' | 'cps450' | 'cps532' | 'cps635'  (see load_input_spot)
    cfg.zemax_unit_to_um   = 1000;    % Zemax POP axes are in lens units (mm -> 1000)
    cfg.comsol_unit_to_um  = [];      % [] = read from COMSOL header / assume um; or 1, 1e-3 (nm)
    cfg.pol_component      = 'Ey';    % input assumed vertically polarized -> Ey ('Ex'|'Ey'|'auto')
    cfg.align_centroids    = true;    % move input centroid onto the mode-1 centroid
    cfg.beam_offset_um     = [0 0];   % extra [dx dy] offset of beam from that position (um)
    cfg.defocus_um         = 0;       % propagate input by this distance (um) before overlap
    cfg.extra_phase_fn     = [];      % optional @(X,Y) radians, X,Y in um relative to beam center
    cfg.lambda_um          = [];      % [] = read from POP header
    cfg.n_medium           = [];      % [] = read from POP header (index of the medium before the facet)
    cfg.phase_sign         = +1;      % flip to -1 if Zemax phase uses the opposite sign convention
    cfg.pad_n              = 2048;    % FFT grid for defocus propagation (zero-padded POP window, ~12 um)
    cfg.grid_n             = 801;     % integration grid points per axis
    cfg.default_dir_irrad  = 'G:\Shared drives\Corp Main\Engineering\Design & SHG Analysis (Brent)\facet_coupling\input_data';
    cfg.default_dir_modes  = 'G:\Shared drives\Corp Main\Engineering\Design & SHG Analysis (Brent)\facet_coupling\input_data';
    cfg.max_montage_modes  = 16;
    % spot_source = 'gaussian': coherent elliptical Gaussian from its mode field diameter
    cfg.gauss = struct('mfd_x_um', 2.0, 'mfd_y_um', 2.0, 'm2', 1, 'lambda_um', 0.45, ...
                       'window_um', 12, 'n_pts', 1024);
    % spot_source = 'measured': intensity image (HDR .mat / tiff / png / csv); see load_measured_spot
    cfg.measured = struct('file', '', 'um_per_pixel', [], 'background', 'auto', 'noise_sigma', 3, ...
                          'crop_um', [6 6], 'rot90', 0, 'flip_x', false, 'flip_y', false, ...
                          'var_name', '', 'lambda_um', 0.45);
    % spot_source = 'cps450' | 'cps532' | 'cps635': through-focus series of that CPS laser's focused
    % beam; the smallest beam is used. (The mode wavelength comes from "wl=<nm>nm" in the mode file name.)
    % (see load_cps_spot). um_per_pixel 0.193 = fit through the origin of 5 calibration targets
    % (individual points scatter by about +-9 %) -- update when the rear-camera calibration is redone.
    % Orientation confirmed (Nov 2026): image vertical = chip y (0.35 um thickness) -> rot90 = 0.
    cfg.cps450 = struct('folder', 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_11_09_LabData\CPS Reference Images\450', ...
                     'um_per_pixel', 0.193, 'z_select', 'auto', 'min_peak', 30, 'sat_level', 250, ...
                     'max_halfwidth_um', 25, 'crop_um', [15 15], 'rot90', 0, 'lambda_um', 0.45, 'plot', true);
    cfg.cps532 = struct('folder', 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_11_09_LabData\CPS Reference Images\532', ...
                     'um_per_pixel', 0.193, 'z_select', 'auto', 'min_peak', 30, 'sat_level', 250, ...
                     'max_halfwidth_um', 25, 'crop_um', [15 15], 'rot90', 0, 'lambda_um', 0.532, 'plot', true);
    cfg.cps635 = struct('folder', 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_11_09_LabData\CPS Reference Images\635', ...
                     'um_per_pixel', 0.193, 'z_select', 'auto', 'min_peak', 30, 'sat_level', 250, ...
                     'max_halfwidth_um', 25, 'crop_um', [15 15], 'rot90', 0, 'lambda_um', 0.635, 'plot', true);
    cfg.mode_lambda_um = [];          % COMSOL mode wavelength (um); used only if the mode file name has no "wl=<nm>nm" tag
    % Optional presets -- leave '' to be prompted with file dialogs
    cfg.irrad_file    = '';           % full path to Zemax POP irradiance .txt
    cfg.phase_file    = '';           % full path to POP phase .txt, or 'none' for flat phase
    cfg.mode_file     = '';           % full path to combined COMSOL mode .csv
    cfg.save_figs_dir = '';           % if set, figures are saved here as PNG

    if nargin > 0 && ~isempty(cfg_override)
        fn = fieldnames(cfg_override);
        for i = 1:numel(fn), cfg.(fn{i}) = cfg_override.(fn{i}); end
    end
end

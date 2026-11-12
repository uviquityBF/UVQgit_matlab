function cfg = modal_default_cfg(cfg_override)
% Default settings shared by modal_power_distribution / modal_power_sweep.
% Any field of the optional struct cfg_override replaces the default.
%
% USAGE:
%   cfg = modal_default_cfg();
%   cfg = modal_default_cfg(struct('mode_file', 'C:\...\modes.csv'));

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

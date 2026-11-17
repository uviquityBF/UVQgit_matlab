function spot = load_input_spot(cfg)
% Loads the input spot from the source chosen in cfg.spot_source and returns one
% struct that load_modal_problem / modal_set_spot work from.
%
%   cfg.spot_source = 'zemax'     Zemax POP irradiance [+ phase] text files (dialogs, or
%                                 cfg.irrad_file / cfg.phase_file)
%                     'gaussian'  coherent Gaussian from its mode field diameter (cfg.gauss)
%                     'measured'  intensity image such as an HDR composite (cfg.measured)
%                     'cps450' | 'cps532' | 'cps635'
%                                 through-focus camera series of that CPS laser; picks the
%                                 smallest beam (cfg.cps450 / cfg.cps532 / cfg.cps635)
%
% Returned fields
%   zX, zY        native-grid axes (um)           E_native  complex field [ny x nx]
%   P_in_total    total power (arb. units)        cx_in, cy_in  intensity centroid (um)
%   lambda_um     wavelength used for defocus     n_medium  refractive index of the medium
%   name          label for titles
% Returns [] if the user cancels a dialog.

    um = char(181);
    src = 'zemax';  if isfield(cfg, 'spot_source') && ~isempty(cfg.spot_source), src = lower(cfg.spot_source); end
    fprintf('--- Step 1: Input spot (%s) ---\n', src);
    switch src
        case 'zemax',    spot = load_zemax_spot(cfg);
        case 'gaussian', spot = make_gaussian_spot(cfg);
        case 'measured', spot = load_measured_spot(cfg);
        case {'cps450', 'cps532', 'cps635'}, spot = load_cps_spot(cfg, src(4:end));
        case 'cps',      spot = load_cps_spot(cfg, '450');      % older name for cps450
        otherwise, error('cfg.spot_source must be ''zemax'', ''gaussian'', ''measured'', ''cps450'', ''cps532'' or ''cps635'' (got "%s").', src);
    end
    if isempty(spot), return; end
    fprintf('Input centroid: (%.2f, %.2f) %sm\n', spot.cx_in, spot.cy_in, um);
    fprintf('Wavelength used for defocus: %.4f %sm in n=%.4f\n', spot.lambda_um, um, spot.n_medium);
end

function spot = load_zemax_spot(cfg)
    start = pwd;  if exist(cfg.default_dir_irrad, 'dir'), start = cfg.default_dir_irrad; end
    if ~isempty(cfg.irrad_file)
        [pathI, nmI, extI] = fileparts(cfg.irrad_file);  fileI = [nmI extI];
    else
        [fileI, pathI] = uigetfile(fullfile(start, '*.txt'), 'Select the Zemax POP IRRADIANCE file');
        if isequal(fileI, 0), disp('Canceled.'); spot = []; return; end
    end
    POP_I = load_zemax_pop_txt(fullfile(pathI, fileI));

    if isempty(cfg.phase_file)
        choice = questdlg('Select a matching Zemax phase file?', 'Phase front', ...
                          'Select phase file', 'Uniform (flat) phase', 'Uniform (flat) phase');
        have_phase = strcmp(choice, 'Select phase file');
    else
        have_phase = ~strcmpi(cfg.phase_file, 'none');
    end
    if have_phase
        if ~isempty(cfg.phase_file)
            [pathP, nmP, extP] = fileparts(cfg.phase_file);  fileP = [nmP extP];
        else
            [fileP, pathP] = uigetfile(fullfile(pathI, '*.txt'), 'Select the matching Zemax POP PHASE file');
            if isequal(fileP, 0), disp('Canceled.'); spot = []; return; end
        end
        POP_P = load_zemax_pop_txt(fullfile(pathP, fileP));
        if ~isequal(size(POP_P.data), size(POP_I.data)) || ...
           abs(POP_P.spacingX - POP_I.spacingX) > 1e-9 * max(1, abs(POP_I.spacingX))
            error('Phase file grid/spacing does not match the irradiance file.');
        end
        phi_native = cfg.phase_sign * POP_P.data;   % radians
        fprintf('Phase: %s\n', fileP);
    else
        phi_native = zeros(size(POP_I.data));
        fprintf('Phase: uniform (flat) -- best-case assumption.\n');
    end

    zX = POP_I.x_axis * cfg.zemax_unit_to_um;      % um
    zY = POP_I.y_axis * cfg.zemax_unit_to_um;
    I_native = max(POP_I.data, 0);
    % Build complex field on the NATIVE grid, so phase wraps are never interpolated.
    E_native = sqrt(I_native) .* exp(1i * phi_native);
    dxz = abs(zX(2) - zX(1));  dyz = abs(zY(2) - zY(1));
    P_in_total = sum(I_native(:)) * dxz * dyz;     % total incident power (arb. units)

    [ZX, ZY] = meshgrid(zX, zY);
    cx_in = sum(ZX(:) .* I_native(:)) / sum(I_native(:));
    cy_in = sum(ZY(:) .* I_native(:)) / sum(I_native(:));

    % Wavelength / medium index for propagation (defocus), from the POP header
    lambda_um = cfg.lambda_um;  n_medium = cfg.n_medium;
    tk = regexp(POP_I.wavelength, 'Wavelength\s+([0-9.eE+-]+)\s+\S+\s+in index\s+([0-9.eE+-]+)', 'tokens', 'once');
    if isempty(lambda_um) && ~isempty(tk), lambda_um = str2double(tk{1}); end
    if isempty(n_medium)  && ~isempty(tk), n_medium  = str2double(tk{2}); end
    if isempty(lambda_um), lambda_um = NaN; end
    if isempty(n_medium),  n_medium = 1; end

    spot.zX = zX;  spot.zY = zY;  spot.E_native = E_native;  spot.P_in_total = P_in_total;
    spot.cx_in = cx_in;  spot.cy_in = cy_in;  spot.lambda_um = lambda_um;  spot.n_medium = n_medium;
    spot.name = fileI;
end

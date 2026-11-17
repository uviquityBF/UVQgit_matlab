function P = load_modal_problem(cfg)
% Loads the input spot (Zemax POP, Gaussian MFD, or measured image; see load_input_spot)
% and a COMSOL mode set,
% places everything on a common integration grid, and returns one struct P that
% modal_power_distribution / modal_power_sweep work from.
%
% USAGE:
%   P = load_modal_problem(modal_default_cfg(...));      % P = [] if canceled
%
% Key fields of P
%   zX, zY, E_native, lambda_um, n_medium, pad_n   input field on its native grid (um)
%   P_in_total, cx_in, cy_in                       input power and centroid
%   intX, intY, IX, IY                             integration grid (um)
%   E_modes, E_other [ny x nx x nM]                mode fields on the grid (input-pol
%                                                  component / the other component)
%   neff, modeIdx, modeTag, polFrac, pol           per-mode info (mode 1 = highest neff)
%   align_shift                                    [dx dy] that puts the input centroid
%                                                  on the mode-1 centroid
%   fileI, modeSetName                             names for titles

    um = char(181);                   % micro sign, avoids file-encoding problems

%% ---------------- 1. LOAD INPUT SPOT ----------------
spot = load_input_spot(cfg);
if isempty(spot), P = []; return; end
zX = spot.zX;  zY = spot.zY;  E_native = spot.E_native;  P_in_total = spot.P_in_total;
cx_in = spot.cx_in;  cy_in = spot.cy_in;  lambda_um = spot.lambda_um;  n_medium = spot.n_medium;
fileI = spot.name;

%% ---------------- 2. LOAD MODE SET ----------------
fprintf('\n--- Step 2: Mode set ---\n');
start = pwd;  if exist(cfg.default_dir_modes, 'dir'), start = cfg.default_dir_modes; end
if ~isempty(cfg.mode_file)
    src = 'Single combined file';
else
    src = questdlg('How are the modes stored?', 'Mode source', ...
                   'Single combined file', 'Directory of files', 'Single combined file');
end
if isempty(src), disp('Canceled.'); P = []; return; end

if strcmp(src, 'Single combined file')
    if ~isempty(cfg.mode_file)
        [pathM, nmM, extM] = fileparts(cfg.mode_file);  fileM = [nmM extM];
    else
        [fileM, pathM] = uigetfile(fullfile(start, '*.csv'), 'Select the combined COMSOL mode export');
        if isequal(fileM, 0), disp('Canceled.'); P = []; return; end
    end
    fprintf('Reading %s ...\n', fileM);
    S = load_comsol_mode_set(fullfile(pathM, fileM), cfg.comsol_unit_to_um);
    nM = numel(S);
    neff = [S.neff]';
    [~, order] = sort(neff, 'descend');           % fundamental (highest neff) first
    S = S(order);  neff = neff(order);
    modeIdx = (1:nM)';
    modeTag = repmat({''}, nM, 1);
    hasEy = [S.hasEy]';
    modes = cell(nM, 1);
    for k = 1:nM, modes{k} = S(k); end
    modeSetName = fileM;
    fprintf('  %d modes, %d points each, neff %.4f ... %.4f\n', nM, numel(S(1).x), neff(1), neff(end));
    mode_lambda_um = parse_mode_wavelength(fileM);
else
    modeDir = uigetdir(start, 'Select directory of COMSOL mode files (one geometry)');
    if isequal(modeDir, 0), disp('Canceled.'); P = []; return; end
    files = [dir(fullfile(modeDir, '*.csv')); dir(fullfile(modeDir, '*.txt'))];
    if isempty(files), error('No .csv/.txt mode files found in %s', modeDir); end

    nM = numel(files);
    modeIdx = nan(nM, 1);  neff = nan(nM, 1);  modeTag = repmat({''}, nM, 1);
    for k = 1:nM
        [~, nm] = fileparts(files(k).name);
        t = regexp(nm, 'mode(\d+)', 'tokens', 'once');
        if ~isempty(t), modeIdx(k) = str2double(t{1}); end
        t = regexp(nm, '_(t[em]\d+)', 'tokens', 'once', 'ignorecase');   % e.g. tm00
        if ~isempty(t), modeTag{k} = upper(t{1}); end
        t = regexp(nm, 'neff(\d+\.?\d*)', 'tokens', 'once');
        if ~isempty(t), neff(k) = str2double(t{1}); end
    end
    % Sort by mode index (files without an index go last, alphabetical)
    sortKey = modeIdx;  sortKey(isnan(sortKey)) = Inf;
    [~, order] = sortrows([sortKey, (1:nM)']);
    files = files(order);  modeIdx = modeIdx(order);  neff = neff(order);  modeTag = modeTag(order);
    for k = 1:nM
        if isnan(modeIdx(k)), modeIdx(k) = k; end
    end

    scale = cfg.comsol_unit_to_um;  if isempty(scale), scale = 1; end
    modes = cell(nM, 1);
    for k = 1:nM
        modes{k} = load_comsol_mode_field(fullfile(modeDir, files(k).name), scale);
        fprintf('  loaded %-32s (%d pts)\n', files(k).name, numel(modes{k}.x));
    end
    hasEy = true(nM, 1);
    modeSetName = modeDir;
    [~, dirName] = fileparts(modeDir);
    mode_lambda_um = parse_mode_wavelength(dirName);
    if isnan(mode_lambda_um), mode_lambda_um = parse_mode_wavelength(files(1).name); end
end

% ---- wavelength the modes were solved at (file-name tag "wl=<nm>nm", else cfg.mode_lambda_um) ----
if isnan(mode_lambda_um)
    mode_lambda_um = getfield_or(cfg, 'mode_lambda_um', NaN);
    if isnan(mode_lambda_um)
        fprintf('  Mode wavelength: unknown (no "wl=<nm>nm" in the file name, cfg.mode_lambda_um empty)\n');
    else
        fprintf('  Mode wavelength: %.0f nm (from cfg.mode_lambda_um)\n', 1000 * mode_lambda_um);
    end
else
    fprintf('  Mode wavelength: %.0f nm (from the file name)\n', 1000 * mode_lambda_um);
end
lam_spot = lambda_um;
if isfield(spot, 'lambda_nominal_um'), lam_spot = spot.lambda_nominal_um; end
if ~isnan(mode_lambda_um) && abs(lam_spot - mode_lambda_um) > 0.01
    warning('load_modal_problem:wavelength', ...
        ['Input spot wavelength (%.0f nm) differs from the wavelength the modes were solved at (%.0f nm): ' ...
         'the overlap is not physically meaningful.'], 1000 * lam_spot, 1000 * mode_lambda_um);
end

%% ---------------- 3. COMMON INTEGRATION GRID ----------------
% Window = extent of the COMSOL cross-section (the chip can only capture
% light inside it; the rest is reported as unaccounted).
xmin = min(modes{1}.x);  xmax = max(modes{1}.x);
ymin = min(modes{1}.y);  ymax = max(modes{1}.y);
intX = linspace(xmin, xmax, cfg.grid_n);
intY = linspace(ymin, ymax, cfg.grid_n);
[IX, IY] = meshgrid(intX, intY);
fprintf('\nIntegration window: X [%.1f, %.1f] %sm, Y [%.1f, %.1f] %sm\n', xmin, xmax, um, ymin, ymax, um);

E_Ex = zeros(cfg.grid_n, cfg.grid_n, nM);
E_Ey = zeros(cfg.grid_n, cfg.grid_n, nM);
F = [];
for k = 1:nM
    m = modes{k};
    sameGrid = ~isempty(F) && isequal(m.x, modes{1}.x) && isequal(m.y, modes{1}.y);
    if sameGrid
        F.Values = real(m.Ex);                  % reuse the triangulation
    else
        F = scatteredInterpolant(m.x, m.y, real(m.Ex), 'linear', 'none');
    end
    re = F(IX, IY);
    F.Values = imag(m.Ex);  im_ = F(IX, IY);
    E_Ex(:, :, k) = nan2zero(re) + 1i * nan2zero(im_);
    F.Values = real(m.Ey);  re = F(IX, IY);
    F.Values = imag(m.Ey);  im_ = F(IX, IY);
    E_Ey(:, :, k) = nan2zero(re) + 1i * nan2zero(im_);
end

% Choose the polarization component used for the (scalar) input field
pExy = [sum(sum(abs(E_Ex(:, :, 1)).^2)), sum(sum(abs(E_Ey(:, :, 1)).^2))];
if strcmpi(cfg.pol_component, 'auto')
    if pExy(1) >= pExy(2), pol = 'Ex'; else, pol = 'Ey'; end
else
    pol = cfg.pol_component;
end
fprintf('Input polarization component matched to mode field: %s\n', pol);
if strcmp(pol, 'Ex'), E_modes = E_Ex; E_other = E_Ey; else, E_modes = E_Ey; E_other = E_Ex; end

% Fraction of each mode's power in the chosen polarization (character check)
polFrac = zeros(nM, 1);
for k = 1:nM
    pa = sum(sum(abs(E_Ex(:, :, k)).^2));  pb = sum(sum(abs(E_Ey(:, :, k)).^2));
    if strcmp(pol, 'Ex'), polFrac(k) = pa / (pa + pb); else, polFrac(k) = pb / (pa + pb); end
end

if strcmp(pol, 'Ey') && any(~hasEy)
    error(['Ey was requested but imag(Ey) is missing from the COMSOL export. ' ...
           'Re-export with imag(ewfd.Ey) as the 4th expression per mode.']);
end

% Auto-tag modes that have no tag: TM if the mode is mostly in the input
% polarization component (Ey -> TM), else TE; then the number of sign changes
% (nodes) of that dominant signed field through the peak along x and y.
% Approximate -- override by naming per-mode files with a _tm00-style tag.
for k = 1:nM
    if isempty(modeTag{k})
        if polFrac(k) >= 0.5, Ed = E_modes(:, :, k); b0 = 'TM'; if strcmp(pol, 'Ex'), b0 = 'TE'; end
        else,                 Ed = E_other(:, :, k); b0 = 'TE'; if strcmp(pol, 'Ex'), b0 = 'TM'; end
        end
        modeTag{k} = sprintf('%s%d%d*', b0, count_nodes(Ed));
    end
end

% ---- alignment of input centroid with the fundamental mode ----
pk = abs(E_modes(:, :, 1)).^2;
cx_m = sum(IX(:) .* pk(:)) / sum(pk(:));
cy_m = sum(IY(:) .* pk(:)) / sum(pk(:));

P = struct('zX', zX, 'zY', zY, 'E_native', E_native, 'lambda_um', lambda_um, ...
    'n_medium', n_medium, 'pad_n', cfg.pad_n, 'P_in_total', P_in_total, ...
    'cx_in', cx_in, 'cy_in', cy_in, 'intX', intX, 'intY', intY, 'IX', IX, 'IY', IY, ...
    'E_modes', E_modes, 'E_other', E_other, 'neff', neff, 'modeIdx', modeIdx, ...
    'modeTag', {modeTag}, 'polFrac', polFrac, 'pol', pol, ...
    'align_shift', [cx_m - cx_in, cy_m - cy_in], 'mode_centroid', [cx_m, cy_m], ...
    'fileI', fileI, 'modeSetName', modeSetName, 'nM', nM, 'mode_lambda_um', mode_lambda_um);
end

function v = getfield_or(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function [nx, ny] = count_nodes(E)
% Number of sign changes of the signed (global-phase-fixed) field along the
% row/column through its peak, ignoring values below 5% of the peak.
    [~, imax] = max(abs(E(:)));
    E = real(E * exp(-1i * angle(E(imax))));
    [r, c] = ind2sub(size(E), imax);
    nx = sign_changes(E(r, :));
    ny = sign_changes(E(:, c));
    if nargout < 2, nx = [nx, ny]; end
end

function n = sign_changes(v)
    v = v(abs(v) > 0.05 * max(abs(v)));
    n = sum(diff(sign(v)) ~= 0);
end

function A = nan2zero(A)
    A(isnan(A)) = 0;
end

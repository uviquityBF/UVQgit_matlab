function [W, gridinfo] = load_mode_width_series(mode_dir, opts)
% Loads a directory of combined COMSOL mode exports, one file per waveguide
% width, onto one common cropped grid ready for eigenmode-expansion overlaps.
%
% File names must contain the width, e.g.  modeData_N=20__h=350nm__w=2500nm.csv
% (pattern  w=<number>nm ). Each file is an Export>Data .csv with 4 columns per
% mode (Re/Im Ex, Re/Im Ey) -- see load_comsol_mode_set. All files must use the
% SAME regular x,y grid (same window and resolution).
%
% USAGE:
%   [W, gridinfo] = load_mode_width_series(dir, opts)
%
% opts fields (all optional)
%   n_cut         modes with Re(neff) below this are non-guided (default 1.8)
%   drop_leaky    true: remove modes below n_cut (their power counts as radiated);
%                 false: keep them and use Im(neff) from the header as loss (default true)
%   crop_x_um     [xmin xmax] crop window (default [-5 5])
%   crop_y_um     [ymin ymax]               (default [-3 3])
%   wavelength_um only load files tagged wl=<nm>nm with this wavelength (untagged files are
%                 kept); required if the folder mixes several wavelengths (default [] = all)
% File names may carry the wavelength the modes were solved at, e.g.
%   modeData_wl=450nm__N=20__h=350nm__w=2500nm.csv   (see parse_mode_wavelength)
%
% W(i), sorted from the WIDEST to the narrowest width:
%   width_um, file
%   A         [2*npts x nk] single complex; each column is one mode, Ex stacked on
%             Ey, normalized so sum(abs(col).^2) = 1 (so overlaps are A2'*A1)
%   neff      [nk x 1] real part;  neff_imag [nk x 1] imaginary part (0 if guided)
%   polFrac   [nk x 1] fraction of mode power in Ey
%   keep_idx  indices of the kept modes within the original descending-neff list
%   n_all, n_dropped, neff_all (Re(neff) of ALL modes in the file, descending)
%   cx_um, cy_um   intensity centroid of the fundamental mode (alignment check)
% gridinfo: struct with x, y (cropped axes, um), nx, ny, lambda_um (mode wavelength from the
%           file names, NaN if untagged)

    if nargin < 2, opts = struct(); end
    n_cut      = getopt(opts, 'n_cut', 1.8);
    drop_leaky = getopt(opts, 'drop_leaky', true);
    crop_x     = getopt(opts, 'crop_x_um', [-5 5]);
    crop_y     = getopt(opts, 'crop_y_um', [-3 3]);

    files = dir(fullfile(mode_dir, '*.csv'));
    width = nan(numel(files), 1);
    for i = 1:numel(files)
        t = regexp(files(i).name, 'w=(\d+\.?\d*)nm', 'tokens', 'once');
        if ~isempty(t), width(i) = str2double(t{1}) / 1000; end
    end
    files = files(~isnan(width));  width = width(~isnan(width));
    if isempty(files), error('No mode files with "w=<number>nm" in the name found in %s', mode_dir); end

    % ---- wavelength tag "wl=<number>nm" (files without a tag are accepted as-is) ----
    wl = arrayfun(@(f) parse_mode_wavelength(f.name), files);
    want = getopt(opts, 'wavelength_um', []);
    if ~isempty(want)
        keep = isnan(wl) | abs(wl - want) < 1e-6;
        if ~any(keep), error('No mode files with wl=%.0fnm in %s', 1000 * want, mode_dir); end
        files = files(keep);  width = width(keep);  wl = wl(keep);
    end
    wl_tags = unique(wl(~isnan(wl)));
    if numel(wl_tags) > 1
        error(['%s holds mode files for several wavelengths (%s nm); set opts.wavelength_um ' ...
               '(cfg.lambda_mode_um in taper_mode_evolution) to pick one.'], mode_dir, num2str(1000 * wl_tags'));
    end
    if isempty(wl_tags), lam_modes = NaN; else, lam_modes = wl_tags; end
    [width, order] = sort(width, 'descend');  files = files(order);
    if any(diff(width) == 0), error('Duplicate widths in %s', mode_dir); end

    W = struct('width_um', cell(numel(files), 1), 'file', [], 'A', [], 'neff', [], ...
               'neff_imag', [], 'polFrac', [], 'keep_idx', [], 'n_all', [], ...
               'n_dropped', [], 'cx_um', [], 'cy_um', [], 'neff_all', []);
    gridinfo = [];
    for i = 1:numel(files)
        fprintf('  [%2d/%2d] w = %.3f um : %s\n', i, numel(files), width(i), files(i).name);
        S = load_comsol_mode_set(fullfile(mode_dir, files(i).name));
        if any(~[S.hasEy]), error('%s lacks imag(Ey); re-export it.', files(i).name); end
        [~, ord] = sort([S.neff], 'descend');  S = S(ord);
        nAll = numel(S);

        % ---- place the scattered list on a regular grid ----
        [xu, ~, ix] = unique(S(1).x);   [yu, ~, iy] = unique(S(1).y);
        if numel(xu) * numel(yu) ~= numel(S(1).x)
            error('%s is not on a regular x,y grid (export with "Regular grid").', files(i).name);
        end
        if isempty(gridinfo)
            cx = xu >= crop_x(1) & xu <= crop_x(2);   cy = yu >= crop_y(1) & yu <= crop_y(2);
            gridinfo.xfull = xu;  gridinfo.yfull = yu;  gridinfo.cx = cx;  gridinfo.cy = cy;
            gridinfo.x = xu(cx);  gridinfo.y = yu(cy);  gridinfo.nx = numel(gridinfo.x);  gridinfo.ny = numel(gridinfo.y);
        elseif numel(xu) ~= numel(gridinfo.xfull) || numel(yu) ~= numel(gridinfo.yfull) || ...
               max(abs(xu - gridinfo.xfull)) > 1e-9 || max(abs(yu - gridinfo.yfull)) > 1e-9
            error('%s uses a different x,y grid than the first file; export all widths on the same x,y grid.', files(i).name);
        end
        lin = sub2ind([numel(yu), numel(xu)], iy, ix);

        npts = gridinfo.nx * gridinfo.ny;
        A = zeros(2 * npts, nAll, 'single');
        for k = 1:nAll
            for c = 1:2
                if c == 1, v = S(k).Ex; else, v = S(k).Ey; end
                G = zeros(numel(yu), numel(xu));  G(lin) = v;
                G = G(gridinfo.cy, gridinfo.cx);
                A((c - 1) * npts + (1:npts), k) = single(G(:));
            end
        end
        nrm = sqrt(sum(abs(A).^2, 1));
        A = bsxfun(@rdivide, A, nrm);

        nre = [S.neff]';  nim = [S.neff_imag]';
        pol = sum(abs(A(npts + 1:end, :)).^2, 1)';          % fraction of power in Ey
        if drop_leaky, keep = find(nre >= n_cut); else, keep = (1:nAll)'; end

        % fundamental-mode centroid (for a consistency check across widths)
        I1 = abs(A(1:npts, 1)).^2 + abs(A(npts + 1:end, 1)).^2;
        [GX, GY] = meshgrid(gridinfo.x, gridinfo.y);
        W(i).cx_um = sum(GX(:) .* double(I1)) / sum(I1);
        W(i).cy_um = sum(GY(:) .* double(I1)) / sum(I1);

        W(i).width_um = width(i);  W(i).file = files(i).name;
        W(i).A = A(:, keep);  W(i).neff = nre(keep);  W(i).neff_imag = nim(keep);
        W(i).polFrac = pol(keep);  W(i).keep_idx = keep;
        W(i).n_all = nAll;  W(i).n_dropped = nAll - numel(keep);  W(i).neff_all = nre;
    end
    gridinfo.lambda_um = lam_modes;       % wavelength from the "wl=<nm>nm" file-name tag (NaN if none)
end

function v = getopt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

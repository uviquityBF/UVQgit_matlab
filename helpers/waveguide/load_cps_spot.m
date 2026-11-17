function spot = load_cps_spot(cfg, which)
% Input spot from a through-focus series of camera images of the focused CPS beam
% (e.g. "CPS Reference Images/450/261111_REF_CPS450+ND40_tint=5000_z=7.94.tiff").
% Every image in the folder is background-subtracted and its second-moment (D4sigma)
% size is measured with an ISO-11146-style iterated integration window. The image with
% the smallest sqrt(D4sigma_x * D4sigma_y) is taken as the waist and handed to
% load_measured_spot (flat phase assumed).
%
% USAGE:  spot = load_cps_spot(cfg, '450')    % which = '450' (default) | '532' | '635'
% The settings come from cfg.cps450 / cfg.cps532 / cfg.cps635 (cfg.cps is accepted as a
% fallback), selected by cfg.spot_source = 'cps450' | 'cps532' | 'cps635'.
%
% Fields of each (see modal_default_cfg)
%   folder        folder with the .tif/.tiff series (one wavelength)
%   um_per_pixel  image scale at the beam plane (required)
%   z_select      'auto' = smallest beam, or a z value (file whose z is closest is used)
%   min_peak      images whose peak is below this many counts above background are ignored
%                 (default 30); images with a peak >= sat_level (default 250) are ignored
%   max_halfwidth_um   cap on the D4sigma integration half-width (default 25) -- keeps the
%                 beam-splitter ghost image and far wings out of the size estimate
%   crop_um       window passed to load_measured_spot for the chosen image (default [15 15])
%   plot          true = through-focus summary figure (default true)
%   rot90, flip_x, flip_y   orientation, passed through to load_measured_spot
% Filenames must contain "tint=<us>" and "z=<number>"; other parts are free.

    if nargin < 2 || isempty(which), which = '450'; end
    which = char(which);
    fld = ['cps' which];
    if isfield(cfg, fld), c = cfg.(fld);
    elseif isfield(cfg, 'cps'), c = cfg.cps;
    else, error('cfg.%s is missing (expected cps450, cps532 or cps635 in modal_default_cfg).', fld);
    end
    if ~isfield(c, 'folder') || isempty(c.folder) || ~exist(c.folder, 'dir')
        error('cfg.%s.folder does not exist: %s', fld, getopt(c, 'folder', '(not set)'));
    end
    % (The spot-vs-mode wavelength check is done in load_modal_problem, which sees the mode file.)
    pix   = getopt(c, 'um_per_pixel', []);
    if isempty(pix), error('cfg.%s.um_per_pixel must be set.', fld); end
    zsel  = getopt(c, 'z_select', 'auto');
    minpk = getopt(c, 'min_peak', 30);
    sat   = getopt(c, 'sat_level', 250);
    hwmax = getopt(c, 'max_halfwidth_um', 25);
    crop  = getopt(c, 'crop_um', [15 15]);
    doplot = getopt(c, 'plot', true);
    um = char(181);

    L = [dir(fullfile(c.folder, '*.tif')); dir(fullfile(c.folder, '*.tiff'))];
    if isempty(L), error('No .tif/.tiff files in %s', c.folder); end
    n = numel(L);
    z = nan(n, 1); tint = nan(n, 1); pk = nan(n, 1); dx = nan(n, 1); dy = nan(n, 1);
    P = nan(n, 1); ok = false(n, 1);
    fprintf('Scanning %d images in %s ...\n', n, c.folder);
    for k = 1:n
        tz = regexp(L(k).name, 'z=([0-9.]+?)\.?tiff?$', 'tokens', 'once');
        tt = regexp(L(k).name, 'tint=([0-9.]+)', 'tokens', 'once');
        if isempty(tz) || isempty(tt), continue; end
        z(k) = str2double(tz{1});  tint(k) = str2double(tt{1});
        A = double(imread(fullfile(c.folder, L(k).name)));
        A = A(:, :, 1);
        [pk(k), dx(k), dy(k), P(k)] = measure_beam(A, pix, hwmax, minpk);
        ok(k) = isfinite(dx(k)) && pk(k) < sat;
    end
    if ~any(ok), error('No usable image (peak between %g counts and saturation).', minpk); end

    d = sqrt(dx .* dy);
    [~, order] = sort(z);
    fprintf('%-62s %7s %7s %9s %9s  %s\n', 'file', 'z', 't(us)', ['D4s_x(' um 'm)'], ['D4s_y(' um 'm)'], 'peak');
    for k = order'
        if isnan(z(k)), continue; end
        fprintf('%-62s %7.3f %7g %9.2f %9.2f  %5.0f%s\n', L(k).name, z(k), tint(k), dx(k), dy(k), pk(k), ...
                ternary(ok(k), '', '  (skipped)'));
    end

    cand = find(ok);
    if ischar(zsel) && strcmpi(zsel, 'auto')
        [~, i] = min(d(cand));
    else
        [~, i] = min(abs(z(cand) - zsel));
    end
    kb = cand(i);
    fprintf('Selected: %s  (z = %.3f, D4sigma %.2f x %.2f %sm)\n', L(kb).name, z(kb), dx(kb), dy(kb), um);
    if ischar(zsel) && strcmpi(zsel, 'auto') && (z(kb) == min(z(cand)) || z(kb) == max(z(cand)))
        warning('load_cps_spot:edge', ...
            'Smallest beam is at the edge of the z range (%.3f): the true waist may lie outside the scan.', z(kb));
    end

    if doplot
        figure('Color', 'w', 'Name', 'CPS through-focus', 'Position', [60 60 760 420]);
        [zs, o] = sort(z(ok));  dxo = dx(ok);  dyo = dy(ok);
        plot(zs, dxo(o), 'o-', 'LineWidth', 1.3); hold on;
        plot(zs, dyo(o), 's-', 'LineWidth', 1.3);
        plot([z(kb) z(kb)], [dx(kb) dy(kb)], 'k*', 'MarkerSize', 10);
        grid on; xlabel('z (stage position, as in file names)'); ylabel(['D4\sigma (' um 'm)']);
        legend({'x', 'y', 'selected'}, 'Location', 'best');
        title(sprintf('CPS beam size through focus (%s, %.3f %sm/px)', L(kb).name, pix, um), 'Interpreter', 'none');
    end

    m = cfg;
    m.measured = struct('file', fullfile(c.folder, L(kb).name), 'um_per_pixel', pix, 'background', 'auto', ...
        'noise_sigma', getopt(c, 'noise_sigma', 3), 'crop_um', crop, 'rot90', getopt(c, 'rot90', 0), ...
        'flip_x', getopt(c, 'flip_x', false), 'flip_y', getopt(c, 'flip_y', false), ...
        'var_name', '', 'lambda_um', getopt(c, 'lambda_um', 0.45));
    spot = load_measured_spot(m);
    spot.name = sprintf('CPS %s', L(kb).name);
end

function [pk, dxum, dyum, ptot] = measure_beam(A, pix, hwmax, minpk)
% Background-subtracted peak and iterated-window D4sigma (um) of one image.
    [nr, nc] = size(A);
    b = max(5, round(0.05 * min(nr, nc)));
    border = [reshape(A(1:b, :), [], 1); reshape(A(end - b + 1:end, :), [], 1); ...
              reshape(A(:, 1:b), [], 1); reshape(A(:, end - b + 1:end), [], 1)];
    bg = median(border);
    sig = 1.4826 * median(abs(border - bg));  if sig <= 0, sig = std(border); end
    pk = max(A(:));                                   % raw peak (saturation check)
    I = A - bg;
    Is = conv2(I, ones(5) / 25, 'same');
    [pkS, ip] = max(Is(:));
    dxum = NaN; dyum = NaN; ptot = NaN;
    if pkS < minpk, return; end
    [r0, c0] = ind2sub(size(I), ip);
    hw = round(hwmax / pix);                          % starting half-width (px)
    for it = 1:8
        rr = max(1, r0 - hw):min(nr, r0 + hw);  cc = max(1, c0 - hw):min(nc, c0 + hw);
        W = I(rr, cc);  W(W < 3 * sig) = 0;
        s = sum(W(:));  if s <= 0, return; end
        [X, Y] = meshgrid(cc, rr);
        cx = sum(X(:) .* W(:)) / s;  cy = sum(Y(:) .* W(:)) / s;
        sx = sqrt(sum((X(:) - cx).^2 .* W(:)) / s);  sy = sqrt(sum((Y(:) - cy).^2 .* W(:)) / s);
        r0 = round(cy);  c0 = round(cx);
        hw = min(round(hwmax / pix), max(10, round(1.5 * 4 * max(sx, sy))));   % 3 x D4sigma window
    end
    dxum = 4 * sx * pix;  dyum = 4 * sy * pix;  ptot = s;
end

function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end

function v = getopt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

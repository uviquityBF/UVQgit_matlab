function spot = load_measured_spot(cfg)
% Imports a measured intensity image of the input spot (e.g. the HDR composite
% from image_analysis/TIFF_HDR_and_Show.m) and converts it to the field form
% used by the modal-decomposition code.
%
% USAGE:
%   spot = load_measured_spot(cfg);     % uses cfg.measured (see below)
%
% Accepted files: .mat (variable I_HDR, or name it in var_name; an optional
% variable um_per_pixel inside the file is used if cfg.measured.um_per_pixel is
% empty), .tif/.tiff/.png/.jpg, .csv (numeric matrix).
%
% cfg.measured fields (all optional)
%   file          path ('' = file dialog)
%   um_per_pixel  image scale at the sample plane. If empty and not in the file, falls
%                 back to 40/360*(50/100) = 0.0556 um/px (the value hard-coded in
%                 input_coupling_overlap.m) with a WARNING: confirm it for your data.
%   background    'auto' = median of a border frame, or a number (counts)
%   noise_sigma   pixels below noise_sigma * (border noise sigma) are set to zero (default 3)
%   crop_um       [half-width_x half-height_y] window around the peak (default [6 6])
%   rot90         number of 90-degree rotations to align image with chip axes (default 0)
%   flip_x, flip_y   mirror the image columns / rows (default false)
%   var_name      variable name inside a .mat ('' = I_HDR, else first 2-D numeric)
%   lambda_um     wavelength, only used if you later defocus (default 0.45)
% Orientation: after rot90/flips, image columns = chip x (lateral), image rows going
% UP = chip +y. The phase is NOT measured: it is assumed flat (best case).
%
% Returns the same struct as load_input_spot.

    m = struct();  if isfield(cfg, 'measured'), m = cfg.measured; end
    file   = getopt(m, 'file', '');
    pix    = getopt(m, 'um_per_pixel', []);
    bgopt  = getopt(m, 'background', 'auto');
    nsig   = getopt(m, 'noise_sigma', 3);
    crop   = getopt(m, 'crop_um', [6 6]);
    nrot   = getopt(m, 'rot90', 0);
    flipx  = getopt(m, 'flip_x', false);
    flipy  = getopt(m, 'flip_y', false);
    vname  = getopt(m, 'var_name', '');
    lam    = getopt(m, 'lambda_um', 0.45);
    um = char(181);

    if isempty(file)
        start = pwd;
        if isfield(cfg, 'default_dir_irrad') && exist(cfg.default_dir_irrad, 'dir'), start = cfg.default_dir_irrad; end
        [f, p] = uigetfile({'*.mat;*.tif;*.tiff;*.png;*.jpg;*.csv', 'Spot image'}, 'Select measured spot image', start);
        if isequal(f, 0), spot = []; return; end
        file = fullfile(p, f);
    end
    [~, nm, ext] = fileparts(file);

    % ---- read ----
    pix_file = [];
    switch lower(ext)
        case '.mat'
            A = load(file);
            if isempty(vname)
                if isfield(A, 'I_HDR'), vname = 'I_HDR';
                else
                    fn = fieldnames(A);  vname = '';
                    for i = 1:numel(fn)
                        v = A.(fn{i});
                        if isnumeric(v) && ismatrix(v) && min(size(v)) > 10, vname = fn{i}; break; end
                    end
                end
            end
            if isempty(vname) || ~isfield(A, vname), error('No 2-D image variable found in %s.', file); end
            I = double(A.(vname));
            if isfield(A, 'um_per_pixel') && isscalar(A.um_per_pixel) && isfinite(A.um_per_pixel) && A.um_per_pixel > 0
                pix_file = A.um_per_pixel;      % NaN / empty in the file means "not set"
            end
        case '.csv'
            I = dlmread(file, ',');
        otherwise
            I = imread(file);
            if size(I, 3) >= 3, I = rgb2gray(I(:, :, 1:3)); end
            I = double(I);
    end
    if isempty(pix), pix = pix_file; end
    if isempty(pix)
        pix = 40 / 360 * (50 / 100);
        warning('load_measured_spot:scale', ...
            ['No um_per_pixel given (cfg.measured.um_per_pixel) or stored in the file; using %.4f um/px ' ...
             'from input_coupling_overlap.m. CONFIRM this calibration for your data.'], pix);
    end

    % ---- orientation ----
    I = rot90(I, nrot);
    if flipx, I = fliplr(I); end
    if flipy, I = flipud(I); end
    [nr, nc] = size(I);

    % ---- background and noise from a border frame ----
    b = max(5, round(0.05 * min(nr, nc)));
    border = [reshape(I(1:b, :), [], 1); reshape(I(end - b + 1:end, :), [], 1); ...
              reshape(I(:, 1:b), [], 1); reshape(I(:, end - b + 1:end), [], 1)];
    if ischar(bgopt), bg = median(border); else, bg = bgopt; end
    sigma = 1.4826 * median(abs(border - median(border)));
    if sigma <= 0, sigma = std(border); end
    I = I - bg;

    % ---- crop around the peak of a lightly smoothed image ----
    Is = conv2(I, ones(3) / 9, 'same');
    [~, ip] = max(Is(:));  [rp, cp] = ind2sub(size(I), ip);
    hx = round(crop(1) / pix);  hy = round(crop(2) / pix);
    rr = rp - hy:rp + hy;  cc = cp - hx:cp + hx;
    C = zeros(numel(rr), numel(cc));
    okr = rr >= 1 & rr <= nr;  okc = cc >= 1 & cc <= nc;
    C(okr, okc) = I(rr(okr), cc(okc));                    % zero-pad if the crop leaves the image
    P_full = sum(max(I(:), 0));

    % ---- remove the noise floor ----
    C(C < nsig * sigma) = 0;
    C = flipud(C);                                        % image rows go down; chip +y goes up
    P_crop = sum(C(:));

    ir = hy + 1;  ic = hx + 1;
    zX = ((1:size(C, 2)) - ic) * pix;
    zY = ((1:size(C, 1)) - ir) * pix;

    spot.zX = zX;  spot.zY = zY;
    spot.E_native = sqrt(C);                              % flat phase
    spot.P_in_total = sum(C(:)) * pix^2;
    [ZX, ZY] = meshgrid(zX, zY);
    spot.cx_in = sum(ZX(:) .* C(:)) / sum(C(:));
    spot.cy_in = sum(ZY(:) .* C(:)) / sum(C(:));
    spot.lambda_um = lam;  spot.n_medium = 1;
    spot.name = [nm ext];

    % ---- report ----
    cx = spot.cx_in;  cy = spot.cy_in;
    sx = sqrt(sum(((ZX(:) - cx).^2) .* C(:)) / sum(C(:)));   % second-moment radii
    sy = sqrt(sum(((ZY(:) - cy).^2) .* C(:)) / sum(C(:)));
    fprintf('Measured spot %s%s: %d x %d px, %.4f %sm/px\n', nm, ext, nc, nr, pix, um);
    fprintf('  background %.2f counts, noise sigma %.2f, peak %.0f counts (SNR %.0f)\n', ...
            bg, sigma, max(Is(:)), max(Is(:)) / max(sigma, eps));
    fprintf('  power kept: %.2f %% of the image total lies inside the %.1f x %.1f %sm crop\n', ...
            100 * P_crop / max(P_full, eps), 2 * crop(1), 2 * crop(2), um);
    fprintf('  FWHM  x = %.3f %sm, y = %.3f %sm\n', hw(sum(C, 1), zX), um, hw(sum(C, 2)', zY), um);
    fprintf('  D4sigma (= MFD for a Gaussian)  x = %.3f %sm, y = %.3f %sm\n', 4 * sx, um, 4 * sy, um);
    fprintf('  Phase is not measured: assuming a flat phase front (best case).\n');
end

function w = hw(v, ax)
% full width at half maximum of a 1-D profile (linear interpolation)
    v = v(:)';  ax = ax(:)';
    [pk, i0] = max(v);  h = pk / 2;
    il = find(v(1:i0) < h, 1, 'last');   ir = i0 - 1 + find(v(i0:end) < h, 1, 'first');
    if isempty(il) || isempty(ir), w = NaN; return; end
    xl = ax(il) + (h - v(il)) * (ax(il + 1) - ax(il)) / (v(il + 1) - v(il));
    xr = ax(ir - 1) + (h - v(ir - 1)) * (ax(ir) - ax(ir - 1)) / (v(ir) - v(ir - 1));
    w = xr - xl;
end

function v = getopt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

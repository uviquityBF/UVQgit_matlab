function Polarization_Camera_Analysis()
% Polarization analysis from a Lucid TRI050S1-P (Sony IMX250MZR-type) micro-
% polarizer camera, exported as four full-resolution JPEGs -- one per polarizer
% orientation (0/45/90/135 deg), already spatially registered by the camera's
% export software.
%
% IMPORTANT PHYSICAL LIMITATION: this sensor only has linear micro-polarizers
% (0/45/90/135 deg). There is no quarter-wave/circular channel, so the circular
% Stokes parameter S3 -- and therefore true ellipticity/handedness -- cannot be
% measured. What this script reports is the LINEAR Stokes vector (S0,S1,S2),
% the degree of LINEAR polarization (DoLP), and the angle of linear polarization
% (AoLP). The "ellipse" plot is the physically correct degenerate case (S3=0):
% a line at angle AoLP, not a true ellipse. If you need ellipticity/handedness,
% you'd need a camera or setup with a rotating quarter-wave plate channel too.
%
% OUTPUTS:
%   (a) a single representative linear-polarization state (S0,S1,S2 -> DoLP,AoLP),
%       averaged in Stokes space over a region of interest (see roi_mode below)
%   (b) per-pixel maps of S0 (total intensity), DoLP, and AoLP -- the "polarization image"
%
% USAGE:
%   1. Set data_dir to the folder with the four *_<angle>.jpg exports.
%   2. Run. Check the raw-intensity figure and the console diagnostics --
%      these images can be very underexposed (dark scenes / heavy ND filtering),
%      in which case most of the frame is noise floor and only a bright spot
%      (e.g. a laser spot) carries a meaningful polarization signal.
%   3. Adjust min_S0_counts / roi_mode as needed and re-run.

    close all; clear; clc;

    %% --- USER SETTINGS ---
    data_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_07_27_LabData\lucid_polarization_camera_test_w_CPS450\';

    min_S0_counts   = 15;    % pixels with S0 below this (8-bit camera counts) are masked as noise-floor / unreliable
    saturation_counts = 250; % warn if any input image has pixels at/above this level (possible clipping)

    roi_mode  = 'auto';      % 'auto' = bright-spot ROI via S0 percentile threshold (good for a laser-spot test)
                              % 'full' = whole frame
                              % 'manual' = use roi_rect below, as [x y w h] in pixels
    roi_percentile = 99.5;   % used when roi_mode = 'auto'
    roi_rect  = [];          % used when roi_mode = 'manual'

    %% --- 1. LOAD AND MATCH THE FOUR ANGLE IMAGES ---
    files = dir(fullfile(data_dir, '*.jpg'));
    if isempty(files)
        error('No .jpg files found in %s', data_dir);
    end

    I0 = []; I45 = []; I90 = []; I135 = [];
    name0 = ''; name45 = ''; name90 = ''; name135 = '';

    fprintf('Loading images from %s ...\n', data_dir);
    for k = 1:numel(files)
        tok = regexp(files(k).name, '_(\d+)[^\d]*\.jpg$', 'tokens', 'once');
        if isempty(tok)
            warning('Skipping %s -- could not find a trailing angle number.', files(k).name);
            continue;
        end
        ang = str2double(tok{1});

        raw = imread(fullfile(data_dir, files(k).name));
        mono = to_mono_image(raw, files(k).name);
        if max(mono(:)) >= saturation_counts
            warning('%s has pixels >= %d counts -- possible clipping.', files(k).name, saturation_counts);
        end

        fprintf('  %-70s -> %g deg channel (max %d counts)\n', files(k).name, ang, max(mono(:)));
        switch ang
            case 0,   I0   = double(mono); name0   = files(k).name;
            case 45,  I45  = double(mono); name45  = files(k).name;
            case 90,  I90  = double(mono); name90  = files(k).name;
            case 135, I135 = double(mono); name135 = files(k).name;
            otherwise
                warning('Skipping %s -- unrecognized angle %g (expected 0/45/90/135).', files(k).name, ang);
        end
    end

    if isempty(I0) || isempty(I45) || isempty(I90) || isempty(I135)
        error('Could not find all four 0/45/90/135 deg images in %s', data_dir);
    end
    if ~isequal(size(I0), size(I45), size(I90), size(I135))
        error('The four images are not all the same size -- check that they came from the same capture.');
    end

    %% --- 2. LINEAR STOKES PARAMETERS ---
    % Averaging both complementary polarizer pairs for S0 halves the shot-noise
    % contribution versus using either pair alone.
    S0 = 0.5 * (I0 + I45 + I90 + I135);
    S1 = I0 - I90;
    S2 = I45 - I135;

    valid = S0 >= min_S0_counts;
    frac_valid = mean(valid(:));
    fprintf('\n%.2f%% of pixels are above the noise floor (S0 >= %g counts).\n', 100*frac_valid, min_S0_counts);
    if frac_valid < 0.02
        warning(['Less than 2%% of the frame has usable signal -- these images look badly underexposed. ' ...
                 'Per-pixel DoLP/AoLP outside the bright region are not meaningful; consider re-shooting ' ...
                 'with more exposure/gain or less ND filtering.']);
    end

    DoLP = nan(size(S0));
    AoLP_deg = nan(size(S0));
    DoLP(valid) = min(1, sqrt(S1(valid).^2 + S2(valid).^2) ./ S0(valid));
    AoLP_deg(valid) = mod(rad2deg(0.5*atan2(S2(valid), S1(valid))), 180); % linear polarization repeats every 180 deg

    %% --- 3. REGION OF INTEREST FOR THE SINGLE REPRESENTATIVE STATE ---
    switch roi_mode
        case 'full'
            mask = true(size(S0));
            roi_label = 'full frame';
        case 'manual'
            if isempty(roi_rect)
                error('roi_mode is "manual" but roi_rect is empty.');
            end
            mask = false(size(S0));
            x = roi_rect(1); y = roi_rect(2); w = roi_rect(3); h = roi_rect(4);
            mask(y:y+h-1, x:x+w-1) = true;
            roi_label = sprintf('manual ROI [x=%d y=%d w=%d h=%d]', roi_rect);
        case 'auto'
            thresh = local_prctile(S0(:), roi_percentile);
            mask = S0 >= thresh;
            roi_label = sprintf('auto ROI (S0 >= %.1fth pct = %.1f counts, %d px)', roi_percentile, thresh, nnz(mask));
        otherwise
            error('Unknown roi_mode "%s"', roi_mode);
    end
    mask = mask & valid;
    if ~any(mask(:))
        error('No valid pixels in the selected ROI -- lower min_S0_counts or roi_percentile.');
    end

    % Average in Stokes space (not per-pixel DoLP/AoLP) -- this is the
    % mathematically correct way to combine partially-polarized measurements.
    S0_roi = mean(S0(mask));
    S1_roi = mean(S1(mask));
    S2_roi = mean(S2(mask));
    DoLP_roi = sqrt(S1_roi^2 + S2_roi^2) / S0_roi;
    AoLP_roi_deg = mod(rad2deg(0.5*atan2(S2_roi, S1_roi)), 180);

    fprintf('\n=== Representative linear polarization state (%s) ===\n', roi_label);
    fprintf('  S0 = %.2f   S1 = %.2f   S2 = %.2f  (counts)\n', S0_roi, S1_roi, S2_roi);
    fprintf('  Degree of linear polarization (DoLP) = %.1f %%\n', 100*DoLP_roi);
    fprintf('  Angle of linear polarization (AoLP)  = %.1f deg (relative to the 0 deg polarizer channel)\n', AoLP_roi_deg);
    fprintf('  NOTE: S3 (circular Stokes parameter) is not measurable with this camera -- \n');
    fprintf('  the ellipse plot below is the degenerate S3=0 case (a line), not a true ellipse.\n\n');

    %% --- 4. PLOTS ---
    plot_intensity_images(I0, I45, I90, I135, name0, name45, name90, name135);
    plot_polarization_maps(S0, DoLP, AoLP_deg, mask);
    plot_representative_ellipse(S0_roi, S1_roi, S2_roi, roi_label);

end


%% ================= Helper functions =================

function mono = to_mono_image(raw, name)
% Extracts a single intensity channel, tolerating RGB-packed grayscale JPEGs
% (common when camera export tools save a mono sensor as a 3-channel JPEG).
    if ismatrix(raw)
        mono = raw;
        return;
    end
    if size(raw,3) ~= 3
        error('%s: unexpected image dimensions.', name);
    end

    r = int16(raw(:,:,1)); g = int16(raw(:,:,2)); b = int16(raw(:,:,3));
    maxdiff = max([max(abs(r-g), [], 'all'), max(abs(r-b), [], 'all')]);
    if maxdiff > 2
        warning(['%s: R/G/B channels differ (max diff %d) -- this looks like a color image, not the ' ...
                  'expected mono-sensor-in-RGB-JPEG. Converting to luma; if this camera is actually a ' ...
                  'color-polarization variant, this script''s polarization math still assumes one ' ...
                  'intensity value per polarizer channel.'], name, maxdiff);
        mono = uint8(0.2989*double(raw(:,:,1)) + 0.5870*double(raw(:,:,2)) + 0.1140*double(raw(:,:,3)));
    else
        mono = raw(:,:,1);
    end
end


function p = local_prctile(x, pct)
% Minimal percentile implementation (avoids requiring the Statistics Toolbox).
    xs = sort(x(:));
    idx = max(1, min(numel(xs), round(pct/100 * numel(xs))));
    p = xs(idx);
end


function plot_intensity_images(I0, I45, I90, I135, name0, name45, name90, name135)
    cmax = max([I0(:); I45(:); I90(:); I135(:)]);
    figure('Name', 'Raw Polarizer-Channel Intensities');
    imgs = {I0, I45, I90, I135};
    labels = {'0 deg', '45 deg', '90 deg', '135 deg'};
    names = {name0, name45, name90, name135};
    for ii = 1:4
        subplot(2,2,ii);
        imagesc(imgs{ii}, [0, cmax]); axis image off; colorbar;
        title({labels{ii}; replace(names{ii}, '_', ' ')}, 'Interpreter', 'none', 'FontSize', 8);
        %colormap(gca, 'gray');
    end
end


function plot_polarization_maps(S0, DoLP, AoLP_deg, mask)
    [ys, xs] = find(mask);
    bbox = [min(xs), min(ys), max(xs)-min(xs), max(ys)-min(ys)];

    figure('Name', 'Polarization Image (per-pixel)');

    subplot(1,3,1);
    imagesc(S0); axis image off; colorbar; colormap(gca, 'gray');
    title('S0: total intensity');
    rectangle('Position', bbox, 'EdgeColor', 'r', 'LineWidth', 1.5);

    subplot(1,3,2);
    im = imagesc(DoLP, [0, 1]); axis image off; colorbar; colormap(gca, 'parula');
    set(im, 'AlphaData', ~isnan(DoLP));
    set(gca, 'Color', [0.5 0.5 0.5]);
    title('DoLP: degree of linear polarization');
    rectangle('Position', bbox, 'EdgeColor', 'r', 'LineWidth', 1.5);

    subplot(1,3,3);
    im = imagesc(AoLP_deg, [0, 180]); axis image off; colorbar; colormap(gca, 'hsv');
    set(im, 'AlphaData', ~isnan(AoLP_deg));
    set(gca, 'Color', [0.5 0.5 0.5]);
    title('AoLP: angle of linear polarization [deg]');
    rectangle('Position', bbox, 'EdgeColor', 'r', 'LineWidth', 1.5);
end


function plot_representative_ellipse(S0_roi, S1_roi, S2_roi, roi_label)
    Ipol = sqrt(S1_roi^2 + S2_roi^2);
    amp_pol = sqrt(max(Ipol, 0));
    amp_tot = sqrt(max(S0_roi, 0));
    aolp = 0.5*atan2(S2_roi, S1_roi);

    theta = linspace(0, 2*pi, 200);
    figure('Name', 'Representative Polarization State');
    hold on;
    plot(amp_tot*cos(theta), amp_tot*sin(theta), 'k--', 'DisplayName', 'total amplitude (S0)');
    plot(amp_pol*[-1, 1]*cos(aolp), amp_pol*[-1, 1]*sin(aolp), 'r-', 'LineWidth', 3, ...
        'DisplayName', 'linear-polarized component');
    axis equal; grid on;
    xlabel('E_x [rel. amplitude]'); ylabel('E_y [rel. amplitude]');
    title({sprintf('Representative state -- %s', roi_label); ...
           sprintf('DoLP = %.1f %%   AoLP = %.1f deg   (S3 not measurable with this sensor)', 100*Ipol/S0_roi, mod(rad2deg(aolp),180))}, ...
           'Interpreter', 'none');
    legend('Location', 'best');
end

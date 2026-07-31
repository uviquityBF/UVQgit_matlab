% Unified TIFF/PNG counting & thresholding tool.
%
% Replaces TIFF_Analysis_Counts_from_Single_Peak.m, TIFF_Analysis_Counts_from_Double_Peak.m,
% and TIFF_Image_Quantification.m (moved to deprecated_scripts\ -- all three did real
% counting/thresholding work but differed only in how the ROI(s) were set up).
%
% For each selected image: flags saturated (>=255) pixels on the raw image, estimates a
% corner-patch background level, sums total image counts, then for each ROI in ROIs applies
% an Otsu threshold and reports the summed/background-subtracted counts and pixel count above
% threshold, with a scatter overlay of the selected pixels.
%
% ROIs        : cell array of ROI rectangles, each [row1 row2 col1 col2]. Pass {} or omit for
%               a single full-frame ROI (reproduces Single_Peak / Quantification behavior).
%               Pass a 2-entry cell array of left/right rectangles to reproduce Double_Peak
%               behavior, e.g. {[400 700 50 500], [400 700 500 950]}.
% thresh_scale: multiplier on the Otsu threshold level (optional, default 1.0)
% start_path  : folder to open the file-picker dialog in (optional)
%
% USAGE:
%   TIFF_Counting_Tool()                                       % single ROI = full frame
%   TIFF_Counting_Tool({[400 700 50 500], [400 700 500 950]})  % two ROIs (left/right)

function TIFF_Counting_Tool(ROIs, thresh_scale, start_path)
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
    close all;

    if nargin < 1 || isempty(ROIs), ROIs = {}; end
    if nargin < 2 || isempty(thresh_scale), thresh_scale = 1.0; end
    if nargin < 3 || isempty(start_path), start_path = pwd; end

    %% Load image stack (kept raw -- saturation is checked before any grayscale conversion)
    [images_raw, filenames, ~] = load_image_stack({'*.tiff';'*.png'}, 'Select image(s)', start_path, 'none');
    Nf = numel(images_raw);

    Nroi = max(1, numel(ROIs));

    %% Loop over files
    for kf = 1:Nf
        I0_raw = images_raw{kf};

        NumSatPix(kf,1) = numel(find(I0_raw >= 255)); %#ok<AGROW>

        I0 = I0_raw;
        if size(I0,3) == 4
            I0 = I0(:,:,1:3);
        end
        if size(I0,3) == 3
            I0 = rgb2gray(I0);
        end

        [BG_avg, BG_std] = estimate_background(I0, [20 20]);
        BGavg_all(kf,1)                 = BG_avg; %#ok<AGROW>
        BG_imagecounts_estimate(kf,1)   = numel(I0) * BG_avg; %#ok<AGROW>
        BG_stdev_estimate(kf,1)         = BG_std; %#ok<AGROW>
        TotalImageCounts(kf,1)          = sum(I0(:)); %#ok<AGROW>
        TOT_stdev_estimate(kf,1)        = std(double(I0(:))); %#ok<AGROW>
        maxval_fullimg                  = max(I0(:));

        Nrows = size(I0,1); Ncols = size(I0,2);

        if isempty(ROIs)
            roi_list = {[1, Nrows, 1, Ncols]};
        else
            roi_list = ROIs;
        end

        hf1 = figure('Name', replace(filenames{kf}, '_', ' ')); %#ok<NASGU>
        roi_str = '';
        for kroi = 1:Nroi
            rect = roi_list{kroi};
            I0roi = I0(rect(1):rect(2), rect(3):rect(4));

            subplot(1, Nroi, kroi);
            [sum_counts, num_pix] = count_roi_above_threshold(I0roi, thresh_scale, true, sprintf('ROI %d', kroi));

            SumCounts(kf,kroi)   = sum_counts; %#ok<AGROW>
            NumPix(kf,kroi)      = num_pix; %#ok<AGROW>
            SumCountsBG(kf,kroi) = sum_counts - num_pix * BG_avg; %#ok<AGROW>

            roi_str = [roi_str, sprintf(' , %d, %d, %d', sum_counts, SumCountsBG(kf,kroi), num_pix)]; %#ok<AGROW>
        end

        fprintf('%s , BGavg=%.2f , maxImg=%d , TotalCounts=%d , NumSatPix=%d%s\n', ...
            filenames{kf}, BG_avg, maxval_fullimg, TotalImageCounts(kf,1), NumSatPix(kf,1), roi_str);
    end

    %% Summary
    fprintf('\n=== Summary (%d file(s), %d ROI(s)) ===\n', Nf, Nroi);
    for kf = 1:Nf
        fprintf('%-40s TotalCounts=%10d  BGavg=%8.2f  NumSatPix=%6d\n', ...
            filenames{kf}, TotalImageCounts(kf,1), BGavg_all(kf,1), NumSatPix(kf,1));
    end
end

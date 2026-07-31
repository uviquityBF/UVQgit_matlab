
% 1. Open a bunch of TIFF files
% 2. Carry out analysis by hand (draw from each section)

function TIFF_Analysis_SHG_from_WG()
%% get data
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
    close all;
    clear;

    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_11_03_LabData\UV_Camera -- SHG from WG';

    %% Load image stack once -- reused by the ratio section and the counting section below
    [I0, filenames, ~] = load_image_stack({'*.tiff';'*.png'}, 'Select image(s)', start_path);
    Nf = numel(I0);

%% Do analyses in sections below ..

 %% [1] Ratio of Two Images
    figure; imagesc(I0{Nf});
    thresh=110;
    isel = find((I0{1}(:) > thresh)& (I0{2}(:) > thresh) & (I0{1}(:)<255));
    isel_not = setdiff( 1:numel(I0{1}), isel);
    Iratio = double(I0{2})./double(I0{1});
    Iratio(isel_not) = 0;
    figure; subplot(1,3,1);  imshow(I0{1}); xlabel(replace(filenames{1},'_',' '));
    subplot(1,3,2);  imshow(I0{2});  xlabel(replace(filenames{2},'_',' '));
    subplot(1,3,3);  imagesc(Iratio); caxis([0.8,1.2]); colorbar;
    [yh,xh] = hist(Iratio(isel),[0.0:0.02:2]);
     figure; bar(xh,yh); xlabel('Ratio of Pixel Value (image 2 /  image 1) = 2xCT450rb / 1xCT450rb'); ylabel('number of pixels')
    title(['median = ',num2str(median(Iratio(isel)))]);

 %% Observe Spacing of other points
 figure;
 for kf=1:min(Nf,6)
    subplot(2,3,kf); imagesc(log10(double(I0{kf})-median(double(I0{kf}))));colorbar
 end


 %% Quantify Bright spot on image
    thresh_scale = 2.0;   % this script's ROIs are the whole frame, thresholded harder than TIFF_Counting_Tool's default
    for kf=1:Nf
        Iimg = I0{kf};

        %-------- TOTAL COUNTS ------------
        [BG_avg, BG_std] = estimate_background(Iimg, [20 20]);
        BG_imagecounts_estimate(kf,1)   = numel(Iimg) * BG_avg; %#ok<AGROW>
        BG_stdev_estimate(kf,1)         = BG_std; %#ok<AGROW>
        TotalImageCounts(kf,1)          = sum(Iimg(:)); %#ok<AGROW>
        TOT_stdev_estimate(kf,1)        = std(double(Iimg(:))); %#ok<AGROW>

        Nrows =  size(Iimg,1);  Ncols = size(Iimg,2);

        %-------- THRESHOLD IMAGE & CALCULATE COUNTS above BG, SHOW FIGURE ------------
        figure;
        [SumCounts(kf,1), NumPix(kf,1)] = count_roi_above_threshold(Iimg, thresh_scale, true, replace(filenames{kf},'_',' ')); %#ok<AGROW>
        [maxval(kf), ~] = max(Iimg(:)); %#ok<AGROW>

        %PRINT RESULT
        X = sprintf('%s , %d , %d, %d, %d, %d, %d, %d, %d, %d',...
            filenames{kf},SumCounts(kf,1),NumPix(kf,1),maxval(kf),...
            TotalImageCounts(kf,1),BG_imagecounts_estimate(kf,1), BG_stdev_estimate(kf,1),TOT_stdev_estimate(kf,1), Nrows, Ncols);
        disp(X);

    end

     disp(maxval)

end

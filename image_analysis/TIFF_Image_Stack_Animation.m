% 1. Open TIFF containing a single peak
%
function TIFF_Image_Stack_Animation()
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
    close all;
    clear;

    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_23_LabData\2024_09_24_UvqA-sd7_Lwgwh520id19.1__SHG_Search_UVImages';

    %% select IMAGE STACK
    [images, filenames, ~] = load_image_stack('*.tiff', 'Select image stack', start_path); %#ok<ASGLU>
    Nf = numel(images);

    %% Read All Files (find center peak of each)
    DAT.I0 = cellfun(@double, images, 'UniformOutput', false);
    for kf = 1:Nf
        DAT.peaks{kf} = find_image_peak(DAT.I0{kf});
    end

    %% select BG IMAGE
    cd(start_path);
    [BG_filename, BG_path] = uigetfile('*.tiff', 'Select background image', 'Multiselect', 'off');
    if isequal(BG_filename, 0)
        error('No background image selected.');
    end

    %% Read BG Image
    BG.I0 = imread(fullfile(BG_path, BG_filename));
    if size(BG.I0,3) == 4
        BG.I0 = BG.I0(:,:,1:3);
    end
    if size(BG.I0,3) == 3
        BG.I0 = rgb2gray(BG.I0);
    end
    BG.I0 = double(BG.I0);
    BG.peak = find_image_peak(BG.I0);

    %% Define ROI
    ROI.r1 = BG.peak.r - 50;
    ROI.r2 = BG.peak.r + 50;
    ROI.c1 = BG.peak.c - 50;
    ROI.c2 = BG.peak.c + 50;

    %% Generate Animation
    hf = figure;
    for kf = 1:Nf
        figure(hf);
        I1 = DAT.I0{kf}(ROI.r1:ROI.r2, ROI.c1:ROI.c2) - BG.I0(ROI.r1:ROI.r2, ROI.c1:ROI.c2);
        imagesc(I1); title(num2str(kf));
        drawnow;
        F(kf) = getframe; %#ok<AGROW>
    end

    figure;
    movie(F,2,1)

end


%% Find center peak of an image via a 10x10 boxcar-smoothed convolution
function peak = find_image_peak(I0)
    A = conv2(I0, ones(10)) / 100;
    [peak.maxval, imax] = max(A(:));
    [peak.r, peak.c] = ind2sub(size(A), imax);
end

% 1. Open stack of "identical" TIFF Images
% 2. Average them
%
function TIFF_Image_Stack_Average()
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
    close all;
    clear;

    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_02_10_LabData\P2P Round2 Measurements';

    %% select IMAGE STACK
    [images, ~, path] = load_image_stack('*.PNG', 'Select image stack', start_path);
    Nf = numel(images);

    %% Read All Files & average
    Iavg = zeros(size(images{1}));
    for kf = 1:Nf
        Iavg = Iavg + double(images{kf});
    end
    Iavg = Iavg / Nf;

    tmp = split(path, '\');
    title_name = replace(tmp{end-1}, '_', ' ');

    figure; imagesc(Iavg); title(title_name);

    save(fullfile(path, 'Iavg.mat'), 'Iavg');
    imwrite(uint16(Iavg), fullfile(path, 'Iavg.tiff'));

    %% Select ROI and analyze
    roi = pick_roi_interactive('Select ROI: click Upper-Left corner', 'Upper-Left selected. Now click Lower-Right corner');
    [IMsel, meanVal, sumVal, Npix] = Analyze_ROI(roi.rect, Iavg, title_name); %#ok<ASGLU>

    %% Shifted ROI (background level, 100 rows down)
    roi_bg = roi;
    roi_bg.rect(1:2) = roi_bg.rect(1:2) + 100;
    [IMsel, meanVal, sumVal, Npix] = Analyze_ROI(roi_bg.rect, Iavg, 'Background Level'); %#ok<ASGLU>

end


%% Analyze ROI
function [IMsel,meanVal,sumVal,Npix] = Analyze_ROI(rect,I0,title_name)
    % compute histogram & average value & STDEV inside ROI
    % rect = [row1 row2 col1 col2]
    IMsel = I0(rect(1):rect(2), rect(3):rect(4));  %select ROI image
    Vsel = IMsel(:);                                    %vector of pixel values
    meanVal = mean(Vsel);
    sumVal = sum(Vsel);
    Npix = length(Vsel);
    [yhist,xhist] = hist(Vsel,40);
    figure;
    subplot(2,1,1);
    imagesc(I0); colorbar;hold on; axis equal;
    plot( [rect(3),rect(4),rect(4),rect(3),rect(3)],...
          [rect(1),rect(1),rect(2),rect(2),rect(1)],'Color','k');
    title({['MeanVal = ',num2str(meanVal)];...
           ['SumVal = ',num2str(sumVal)];...
           ['NPIX = ',num2str(Npix)];...
            title_name});
    subplot(2,1,2); bar(xhist,yhist); grid on;

    disp([num2str(meanVal),',',num2str(sumVal),',',num2str(Npix)]);
end

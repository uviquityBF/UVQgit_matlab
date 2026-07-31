% Open Image, Select ROI,  Deliver MONOCHROME Statistics
% 18 July 2024

%%
addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
clear;

program_name = 'Get Image Statistics';

%select file
dirs.d = 'C:\Users\brent\MATLAB\matDat';
cd(dirs.d); [FileName,PathName] = uigetfile({'*.tiff'},program_name);
filenametitle = replace(FileName,'_',' ');

%open file
if strcmp(FileName(end-3:end),'tiff')  %for TIFF file input
    I0 = imread([PathName,FileName]);
    [Nr,Nc,Nlayers] = size(I0);
    if size(I0,3)>=3
        I1 = rgb2gray(I0); %convert to monochrom
    else
        I1 = I0;        %take as monochrome
    end
elseif strcmp(FileName(end-3:end),'.mat')   %for HDR (*.mat) file input
    INDAT = load([PathName,FileName]);
    I0 = uint32( INDAT.I_HDR );
    [Nr,Nc,Nlayers] = size(I0);
    I1 = I0(:,:,1);
end

%show monochrome image and define ROI corners
hf1 = figure;
imagesc(log10(double(I1))); colorbar;
roi = pick_roi_interactive('Click first corner of ROI', 'Click opposite corner of ROI');

%Get ROI and compute Statistics
I_roi = I1(roi.y(1):roi.y(2), roi.x(1):roi.x(2));

I_roi_mean = mean(I_roi(:));
I_roi_stdev = std(double(I_roi(:)));
disp(FileName);
disp(I_roi_mean);
disp(I_roi_stdev);

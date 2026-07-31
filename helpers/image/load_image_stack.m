function [images, filenames, folder] = load_image_stack(filter_spec, dialog_title, start_path, convert_mode)
% Interactive multi-file picker + loader for a stack of TIFF/PNG images.
%
% filter_spec : uigetfile filter spec, e.g. {'*.tiff';'*.png'} or '*.tiff'
% dialog_title: dialog box title (optional, default '')
% start_path  : folder to cd into before opening the dialog (optional, default pwd)
% convert_mode: 'gray' (default) drops any alpha channel and converts a
%               3-channel image to grayscale; 'none' returns each image
%               exactly as imread loaded it (needed when a caller wants to
%               inspect individual color channels, e.g. a Blue-channel ROI).
%
% Returns:
%   images    : 1xN cell array of images, native class as loaded
%   filenames : 1xN cellstr of selected filenames (no path)
%   folder    : folder the files were selected from, as returned by uigetfile
%
% Shared by TIFF_Analysis_Counts_from_Single_Peak.m, TIFF_Analysis_Counts_from_Double_Peak.m,
% TIFF_Analysis_SHG_from_WG.m, TIFF_Image_Stack_Animation.m, TIFF_Image_Stack_Average.m, and
% define_ROI_for_P2P_image.m, which each independently implemented this uigetfile-Multiselect
% + cellstr-normalize + imread loop.

    if nargin < 2 || isempty(dialog_title), dialog_title = ''; end
    if nargin < 3 || isempty(start_path), start_path = pwd; end
    if nargin < 4 || isempty(convert_mode), convert_mode = 'gray'; end

    if exist(start_path, 'dir')
        cd(start_path);
    end

    [filenames, folder] = uigetfile(filter_spec, dialog_title, 'Multiselect', 'on');
    if isequal(filenames, 0)
        error('load_image_stack:noSelection', 'No file(s) selected.');
    end
    if ischar(filenames)
        filenames = {filenames};
    end
    Nf = numel(filenames);

    images = cell(1, Nf);
    for kf = 1:Nf
        I0 = imread(fullfile(folder, filenames{kf}));
        if strcmp(convert_mode, 'gray')
            if size(I0,3) == 4
                I0 = I0(:,:,1:3);
            end
            if size(I0,3) == 3
                I0 = rgb2gray(I0);
            end
        end
        images{kf} = I0;
    end

    fprintf('Loaded %d image(s) from %s\n', Nf, folder);
end

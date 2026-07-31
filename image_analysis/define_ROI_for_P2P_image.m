
function A = define_ROI_for_P2P_image( )
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\image');
    %clear & get files
    clear; close all;
    program_name = 'define ROI for P2P';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL00\am123\sd01\Standardized Testing';

    [images, filenames, PathName] = load_image_stack({'*.tiff'}, program_name, dirs.d, 'none');
    Nf = numel(images);

    fID = fopen(fullfile(PathName, 'tmp_ROIs.csv'), 'w');
    for kf = 1:Nf
        output_str = select_ROI_from_Image(images{kf}, filenames{kf});
        fprintf(fID,[output_str,'\r\n']);
    end
    fclose(fID);

end

function output_str = select_ROI_from_Image(I0, FileName)
     %take Blue layer of a color TIFF, or the image as-is if monochrome
    if size(I0,3) >= 3
        I1 = I0(:,:,3);
    else
        I1 = I0;
    end

    % show image and select ROI
    hf = figure; imagesc(log10(double(I1)));
    roi = pick_roi_interactive('Click opposite Corners of ROI', 'Now click the other corner');

    output_str = sprintf('%s, %d, %d, %d, %d', FileName, roi.x(1), roi.y(1), roi.x(2), roi.y(2));
    disp(output_str)
    close(hf);
end

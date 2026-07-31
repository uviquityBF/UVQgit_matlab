function roi = pick_roi_interactive(msg1, msg2)
% Prompts the user to click two opposite corners of a rectangular ROI on
% the CURRENT figure/axes (display the image first), draws the resulting
% rectangle, and returns sorted, rounded pixel coordinates.
%
% msg1, msg2: message-box text shown before the first/second click (optional)
%
% Returns roi, a struct with:
%   .x, .y : [min max] sorted, rounded pixel coordinates of the two clicks
%   .rect  : [row1 row2 col1 col2], ready for I(row1:row2, col1:col2) indexing
%
% Shared by Get_Image_Statistics.m, TIFF_Image_Stack_Average.m, and
% define_ROI_for_P2P_image.m, which each independently implemented a
% msgbox + two ginput calls + corner-sort ROI picker.

    if nargin < 1 || isempty(msg1), msg1 = 'Click first corner of ROI'; end
    if nargin < 2 || isempty(msg2), msg2 = 'Click opposite corner of ROI'; end

    h1 = msgbox(msg1);
    uiwait(h1);
    [x1, y1] = ginput(1);
    hold on;
    plot(x1, y1, '+', 'Color', 'y');

    h2 = msgbox(msg2);
    uiwait(h2);
    [x2, y2] = ginput(1);
    plot(x2, y2, '+', 'Color', 'y');
    plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 'Color', 'y');

    roi.x = round(sort([x1, x2]));
    roi.y = round(sort([y1, y2]));
    roi.rect = [roi.y(1), roi.y(2), roi.x(1), roi.x(2)];
end

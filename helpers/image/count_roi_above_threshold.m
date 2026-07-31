function [sum_counts, num_pix, sel_mask, level] = count_roi_above_threshold(I_roi, thresh_scale, show_plot, plot_title)
% Otsu-thresholds an ROI image and sums/counts the pixels above the scaled
% threshold. Optionally overlays the selected pixels as a red scatter on
% top of the ROI image (on the CURRENT axes -- set up figure/subplot first).
%
% I_roi       : 2-D ROI image, any numeric class
% thresh_scale: multiplier on the Otsu (graythresh) level (optional, default 1.0)
% show_plot   : true/false, overlay selected pixels on the current axes (optional, default true)
% plot_title  : title applied to the plot when show_plot is true (optional)
%
% Returns:
%   sum_counts: sum of I_roi pixel values above threshold
%   num_pix   : number of pixels above threshold
%   sel_mask  : logical mask of selected (above-threshold) pixels, same size as I_roi
%   level     : the raw (unscaled) Otsu threshold level used
%
% Shared by TIFF_Analysis_Counts_from_Single_Peak.m, TIFF_Analysis_Counts_from_Double_Peak.m,
% and TIFF_Analysis_SHG_from_WG.m, which each independently ran graythresh/im2bw on an ROI and
% overlaid the selected pixels.

    if nargin < 2 || isempty(thresh_scale), thresh_scale = 1.0; end
    if nargin < 3 || isempty(show_plot), show_plot = true; end
    if nargin < 4, plot_title = ''; end

    level = graythresh(I_roi);
    sel_mask = im2bw(I_roi, level * thresh_scale); %#ok<IM2BW> -- kept for R2018b compatibility

    sum_counts = sum(I_roi(sel_mask));
    num_pix = sum(sel_mask(:));

    if show_plot
        [rsel, csel] = find(sel_mask);
        imshow(I_roi); hold on;
        plot(csel, rsel, '.', 'Color', 'r');
        if ~isempty(plot_title)
            title(plot_title);
        end
    end
end

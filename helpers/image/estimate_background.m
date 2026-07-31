function [bg_mean, bg_std, patch] = estimate_background(I0, patch_size)
% Estimates background level from a corner patch (top-left) of an image.
%
% I0        : 2-D image, any numeric class
% patch_size: [rows, cols] corner patch size (optional, default [20 20]),
%             clipped to the image dimensions if the image is smaller
%
% Returns:
%   bg_mean: mean pixel value in the corner patch (double)
%   bg_std : standard deviation of pixel values in the corner patch (double)
%   patch  : the corner patch itself, for further inspection if desired
%
% Shared by TIFF_Analysis_Counts_from_Single_Peak.m, TIFF_Analysis_Counts_from_Double_Peak.m,
% TIFF_Analysis_SHG_from_WG.m, and TIFF_Image_Quantification.m, which each independently took a
% small corner patch (20x20, or 10x10 in TIFF_Image_Quantification) as a background estimate.

    if nargin < 2 || isempty(patch_size), patch_size = [20, 20]; end
    if isscalar(patch_size), patch_size = [patch_size, patch_size]; end

    nr = min(patch_size(1), size(I0,1));
    nc = min(patch_size(2), size(I0,2));
    patch = I0(1:nr, 1:nc);

    bg_mean = mean(double(patch(:)));
    bg_std  = std(double(patch(:)));
end

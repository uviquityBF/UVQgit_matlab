function [center_wl, maxVal, avg_noise, wl_seg, int_seg, is_valid] = spectral_centroid(wl_raw, int_raw, wl_min, wl_max, varargin)
% Adaptive-threshold spectral centroid (peak wavelength) from one raw spectrum.
%
% Baseline noise is estimated over baseline_band_nm and used to build an
% adaptive detection threshold (avg_noise + noise_k*std_noise). If the peak
% within [wl_min, wl_max] clears that threshold and isn't sitting at the edge
% of the search window, center_wl is the intensity-weighted mean wavelength
% of samples above a cutoff between the baseline and the peak (centroid_frac).
% Otherwise center_wl = NaN and is_valid = false.
%
% wl_seg/int_seg (the raw segment masked to [wl_min, wl_max]) are returned so
% callers can build the same normalized-overlay plot both original scripts did,
% without re-deriving the mask.
%
% Shared by Analyze_LIV_Spectra.m and Analyze_LIV_and_Spectra_Final_PowerRamp.m,
% which independently implemented this identical logic.
%
% Optional name-value pairs (defaults match both original callers):
%   'BaselineBandNm'        [lo hi] nm used to estimate avg_noise/std_noise (default [200 300])
%   'NoiseK'                auto_thresh = avg_noise + NoiseK*std_noise (default 10)
%   'CentroidFrac'          cutOff = avg_noise + CentroidFrac*(maxVal-avg_noise) (default 0.5)
%   'EmptyBaselineFallback' [avg_noise, auto_thresh] used when BaselineBandNm contains
%       no samples of wl_raw (default [0, 100]). NOTE: Analyze_LIV_Spectra.m had this
%       guard; Analyze_LIV_and_Spectra_Final_PowerRamp.m did not (it would have hit an
%       error computing mean([]) on an empty baseline band). Folding both callers onto
%       this shared guard is a (harmless) behavior change for that script -- flagged here.

    p = inputParser;
    p.addParameter('BaselineBandNm', [200, 300]);
    p.addParameter('NoiseK', 10);
    p.addParameter('CentroidFrac', 0.5);
    p.addParameter('EmptyBaselineFallback', [0, 100]);
    p.parse(varargin{:});
    opt = p.Results;

    center_wl = NaN;
    maxVal = NaN;
    wl_seg = [];
    int_seg = [];
    is_valid = false;

    baseMask = (wl_raw >= opt.BaselineBandNm(1)) & (wl_raw <= opt.BaselineBandNm(2));
    if any(baseMask)
        avg_noise = mean(int_raw(baseMask));
        std_noise = std(int_raw(baseMask));
        auto_thresh = avg_noise + opt.NoiseK * std_noise;
    else
        avg_noise = opt.EmptyBaselineFallback(1);
        auto_thresh = opt.EmptyBaselineFallback(2);
    end

    mask = (wl_raw >= wl_min) & (wl_raw <= wl_max);
    wl_seg = wl_raw(mask);
    int_seg = int_raw(mask);
    if isempty(wl_seg)
        return;
    end

    [maxVal, maxIdx] = max(int_seg);
    if (maxVal > auto_thresh) && (maxIdx > 1) && (maxIdx < length(int_seg))
        cutOff = avg_noise + (maxVal - avg_noise) * opt.CentroidFrac;
        peakMask = int_seg >= cutOff;
        center_wl = sum(wl_seg(peakMask) .* int_seg(peakMask)) / sum(int_seg(peakMask));
        is_valid = true;
    end
end

function fwhm_dk = shg_compute_pm_fwhm(dk_vec, shg_norm)
% Full-width-half-max of a normalized phase-matching tuning curve.
    fwhm_dk = NaN;
    [~, peak_idx] = max(shg_norm);
    half_max = 0.5;

    t_left = NaN;
    for i = peak_idx:-1:2
        if shg_norm(i-1) <= half_max && shg_norm(i) >= half_max
            t_left = interp1([shg_norm(i-1), shg_norm(i)], [dk_vec(i-1), dk_vec(i)], half_max);
            break;
        end
    end

    t_right = NaN;
    for i = peak_idx:(length(shg_norm)-1)
        if shg_norm(i) >= half_max && shg_norm(i+1) <= half_max
            t_right = interp1([shg_norm(i), shg_norm(i+1)], [dk_vec(i), dk_vec(i+1)], half_max);
            break;
        end
    end

    if ~isnan(t_left) && ~isnan(t_right)
        fwhm_dk = t_right - t_left;
    end
end

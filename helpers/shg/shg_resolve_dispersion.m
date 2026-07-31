function [dn_pump, dn_shg] = shg_resolve_dispersion(ref_case, json_path)
% Resolves dn/dlambda for pump and SH from a COMSOL dispersion-table JSON
% (interpolated by geometry), falling back to hardcoded median values with an
% interactive override prompt if the table is missing or unreadable.
    MEDIAN_PUMP = -0.001476;
    MEDIAN_SHG  = -0.013490;

    h_nm = ref_case.height_um * 1e3;
    w_nm = ref_case.width_um  * 1e3;

    if exist(json_path, 'file')
        try
            raw  = jsondecode(fileread(json_path));
            data = raw.data;
            h_t  = double([data.h_core_nm]');
            w_t  = double([data.w_core_nm]');
            dp_t = double([data.dn_dWL_fund]');
            ds_t = double([data.dn_dWL_target]');

            Fp = scatteredInterpolant(h_t, w_t, dp_t, 'linear', 'linear');
            Fs = scatteredInterpolant(h_t, w_t, ds_t, 'linear', 'linear');
            dn_pump = Fp(h_nm, w_nm);
            dn_shg  = Fs(h_nm, w_nm);

            in_bounds = shg_check_table_bounds(h_t, w_t, h_nm, w_nm);
            fprintf('\nDispersion resolved from JSON  (h=%.0fnm, w=%.0fnm):\n', h_nm, w_nm);
            fprintf('  dn_pump/d\x03BB = %.6f /nm\n', dn_pump);
            fprintf('  dn_SH/d\x03BB   = %.6f /nm\n', dn_shg);
            if ~in_bounds
                fprintf('  WARNING: geometry outside table range (h=335-365nm, w=350-550nm) — extrapolated.\n');
            end
            return;
        catch err
            fprintf('\nWARNING: dispersion_table.json found but could not be read: %s\n', err.message);
        end
    else
        fprintf('\ndispersion_table.json not found in script directory.\n');
    end

    fprintf('Hardcoded median values (COMSOL table, 21 geometries):\n');
    fprintf('  dn_pump/d\x03BB = %.6f /nm\n', MEDIAN_PUMP);
    fprintf('  dn_SH/d\x03BB   = %.6f /nm\n', MEDIAN_SHG);

    resp = input('Accept these values? [Y/n]: ', 's');
    if strcmpi(strtrim(resp), 'n') || strcmpi(strtrim(resp), 'no')
        val = input(sprintf('  Enter dn_pump/d\x03BB [1/nm, default %.6f]: ', MEDIAN_PUMP), 's');
        dn_pump = str2double(strtrim(val));
        if isnan(dn_pump), dn_pump = MEDIAN_PUMP; end

        val = input(sprintf('  Enter dn_SH/d\x03BB   [1/nm, default %.6f]: ', MEDIAN_SHG), 's');
        dn_shg = str2double(strtrim(val));
        if isnan(dn_shg), dn_shg = MEDIAN_SHG; end
    else
        dn_pump = MEDIAN_PUMP;
        dn_shg  = MEDIAN_SHG;
    end

    fprintf('Using: dn_pump/d\x03BB = %.6f /nm,  dn_SH/d\x03BB = %.6f /nm\n', dn_pump, dn_shg);
end


function in_bounds = shg_check_table_bounds(h_t, w_t, h, w)
    in_bounds = (h >= min(h_t)) && (h <= max(h_t)) && ...
                (w >= min(w_t)) && (w <= max(w_t));
end

function cases = shg_resolve_taper_params(cases, dk_tune)
% RESOLVE TAPER PARAMETERS
%
% Single-segment (no L1_um in taper struct):
%   dlam_PM_dz = (dlam_PM_dh*dh_per_mm + dlam_PM_dw*dw_per_mm) * 1e-3  [nm/um]
%   sets cases(k).dk_center and cases(k).dk_grad
%
% Two-segment, width-endpoint spec (preferred): L1_um + w_start_nm/w_L1_nm/w_end_nm
%   Seg1 rate derived from (w_L1_nm - w_start_nm) / L1_um -- fixed since L1 is constant.
%   Seg2 rate deferred: stored as dk_grad_h + C_scale_dw factors so shg_run_core_simulation
%   can recompute dk_grad_2 per-L (seg2 length = L - L1 varies across the sweep).
%   Sets: cases(k).dk_center, dk_grad (seg1), L1_um, w_L1_nm, w_end_nm, dk_grad_h, C_scale_dw.
%
% Two-segment, rate spec (legacy): L1_um + dw_per_mm_1 + dw_per_mm_2
%   Both dk_grad values pre-computed; sets cases(k).dk_grad_2 directly.
    any_taper = false;
    for k = 1:length(cases)
        if ~isfield(cases(k), 'taper') || isempty(cases(k).taper)
            continue;
        end
        t = cases(k).taper;
        if ~isfield(t, 'lam_PM_0_nm') || isnan(t.lam_PM_0_nm)
            continue;
        end
        if ~any_taper
            fprintf('\nResolving geometry-based chirp parameters...\n');
            any_taper = true;
        end

        lam_nm    = cases(k).lam_nm;
        C_rad_nm2 = (4*pi/lam_nm) * (dk_tune.dn_shg_per_nm/2 - dk_tune.dn_pump_per_nm);
        C_scale   = C_rad_nm2 * 1000;   % rad/(um*nm)

        dk_ctr = C_scale * (lam_nm - t.lam_PM_0_nm);

        if isfield(t, 'L1_um') && ~isempty(t.L1_um)
            % --- Two-segment width taper ---
            dlam_dz_h = t.dlam_PM_dh * t.dh_per_mm * 1e-3;   % height contribution [nm/um], global

            if isfield(t, 'w_start_nm')
                % Width-endpoint spec: actual widths at input, boundary, and output.
                % Seg1 rate is fixed (L1 is constant); seg2 rate is deferred to shg_run_core_simulation.
                dw_mm_1   = (t.w_L1_nm - t.w_start_nm) / (t.L1_um * 1e-3);   % [nm/mm]
                dlam_dz_1 = dlam_dz_h + t.dlam_PM_dw * dw_mm_1 * 1e-3;
                dk_grd_1  = -C_scale * dlam_dz_1;

                cases(k).dk_center  = dk_ctr;
                cases(k).dk_grad    = dk_grd_1;
                cases(k).L1_um      = t.L1_um;
                cases(k).w_L1_nm    = t.w_L1_nm;
                cases(k).w_end_nm   = t.w_end_nm;
                cases(k).dk_grad_h  = -C_scale * dlam_dz_h;         % height-only dk gradient [rad/um^2]
                cases(k).C_scale_dw = -C_scale * t.dlam_PM_dw * 1e-3;  % width factor [rad/um^2 per nm/mm]

                fprintf('  Case %d (%s):  [two-segment, width endpoint]\n', k, cases(k).name);
                fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
                fprintf('    Seg1 (L1=%d um): w %.0f->%.0f nm  dw=%.2f nm/mm  =>  dk_grad_1 = %.4e rad/um^2\n', ...
                        t.L1_um, t.w_start_nm, t.w_L1_nm, dw_mm_1, dk_grd_1);
                fprintf('    Seg2: w %.0f->%.0f nm  (dk_grad_2 computed per L at runtime)\n', ...
                        t.w_L1_nm, t.w_end_nm);
                fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
            else
                % Legacy rate spec: dw_per_mm_1 and dw_per_mm_2 (pre-computable).
                dlam_dz_1 = dlam_dz_h + t.dlam_PM_dw * t.dw_per_mm_1 * 1e-3;
                dlam_dz_2 = dlam_dz_h + t.dlam_PM_dw * t.dw_per_mm_2 * 1e-3;
                dk_grd_1  = -C_scale * dlam_dz_1;
                dk_grd_2  = -C_scale * dlam_dz_2;

                cases(k).dk_center = dk_ctr;
                cases(k).dk_grad   = dk_grd_1;
                cases(k).dk_grad_2 = dk_grd_2;
                cases(k).L1_um     = t.L1_um;

                fprintf('  Case %d (%s):  [two-segment, rate spec]\n', k, cases(k).name);
                fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
                fprintf('    Seg1 (L1=%d um):   dlam_PM/dz = %.4f nm/mm  =>  dk_grad_1 = %.4e rad/um^2\n', t.L1_um, dlam_dz_1*1e3, dk_grd_1);
                fprintf('    Seg2 (remainder): dlam_PM/dz = %.4f nm/mm  =>  dk_grad_2 = %.4e rad/um^2\n', dlam_dz_2*1e3, dk_grd_2);
                fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
            end
        else
            % --- Single segment (v3.2 behavior) ---
            dlam_PM_dz = (t.dlam_PM_dh * t.dh_per_mm + t.dlam_PM_dw * t.dw_per_mm) * 1e-3;
            dk_grd     = -C_scale * dlam_PM_dz;

            cases(k).dk_center = dk_ctr;
            cases(k).dk_grad   = dk_grd;

            fprintf('  Case %d (%s):\n', k, cases(k).name);
            fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
            fprintf('    dlam_PM/dz = %.4f nm/mm  =>  dk_grad = %.4e rad/um^2\n', dlam_PM_dz*1e3, dk_grd);
            fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
        end
    end
end

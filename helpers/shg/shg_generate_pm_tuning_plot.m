function shg_generate_pm_tuning_plot(cases, dk_tune, sim, blue_losses)
% Phase-matching tuning curves (normalized SHG efficiency vs. pump-wavelength
% detuning) for a reference case, at several waveguide lengths and loss levels.
    ref_case = cases(end);

    lam_nm    = ref_case.lam_nm;
    dn_pump   = dk_tune.dn_pump_per_nm;
    dn_shg    = dk_tune.dn_shg_per_nm;
    C_rad_nm2 = (4*pi/lam_nm) * (dn_shg/2 - dn_pump);
    C_scale   = C_rad_nm2 * 1000;

    [a0_sorted, ~] = sort(blue_losses(dk_tune.loss_idx), 'ascend');
    n_loss_sel     = length(a0_sorted);

    cfg.sim              = sim;
    cfg.sim.dz           = dk_tune.dz_um;
    cfg.phys             = ref_case;
    cfg.phys.Pp_peak_mW  = cfg.phys.Pp_avg_mW / cfg.phys.duty_factor;
    cfg.phys.g           = shg_g_coefficient(cfg);
    cfg.phys.a3_abs_dBcm = 0;

    L_vals  = dk_tune.L_vals_um;
    n_L     = length(L_vals);
    n_pts   = dk_tune.n_pts;
    colors  = lines(n_L);
    dk_max  = abs(C_scale) * dk_tune.dl_max_nm;
    dk_vec  = linspace(-dk_max, dk_max, n_pts);

    all_dl_nm  = cell(n_L, n_loss_sel);
    all_guided = cell(n_L, n_loss_sel);
    all_scat   = cell(n_L, n_loss_sel);
    fwhm_g_nm  = NaN(n_L, n_loss_sel);
    fwhm_s_nm  = NaN(n_L, n_loss_sel);

    for lj = 1:n_loss_sel
        a0_ref = a0_sorted(lj);
        cfg.phys.a0_dBcm      = a0_ref;
        cfg.phys.a3_scat_dBcm = a0_ref * ref_case.uv_loss_val;
        fprintf('\n  Loss level: a0=%d dB/cm\n', a0_ref);

        for li = 1:n_L
            L = L_vals(li);
            guided_raw = zeros(1, n_pts);
            scat_raw   = zeros(1, n_pts);
            fprintf('    L=%d um (%d points)...\n', L, n_pts);

            for di = 1:n_pts
                cfg.phys.L_um      = L;
                cfg.phys.dk_center = dk_vec(di);
                [Pg, Ps] = shg_run_core_simulation(cfg);
                guided_raw(di) = Pg * ref_case.duty_factor * ref_case.OCE;
                scat_raw(di)   = Ps * ref_case.duty_factor * ref_case.OCE;
            end

            pk_g  = max(guided_raw);
            pk_s  = max(scat_raw);
            dl_nm = dk_vec / C_scale;

            all_dl_nm{li, lj}  = dl_nm;
            all_guided{li, lj} = guided_raw / pk_g;
            all_scat{li, lj}   = scat_raw   / pk_s;
            fwhm_g_nm(li, lj)  = shg_compute_pm_fwhm(dl_nm, all_guided{li, lj});
            fwhm_s_nm(li, lj)  = shg_compute_pm_fwhm(dl_nm, all_scat{li, lj});
        end
    end

    disp_note = sprintf('dn_{pump}/d\\lambda=%.5f /nm,  dn_{SH}/d\\lambda=%.5f /nm', dn_pump, dn_shg);

    figure('Color','w','Name','Phase Matching Tuning Curves', ...
           'Position',[250, 60, 920, 320*n_loss_sel]);

    for lj = 1:n_loss_sel
        ax = subplot(n_loss_sel, 1, lj);
        hold(ax, 'on');

        h_leg      = gobjects(n_L, 1);
        leg_labels = cell(n_L, 1);

        for li = 1:n_L
            h_leg(li) = plot(ax, all_dl_nm{li,lj}, all_guided{li,lj}, '-', ...
                             'Color', colors(li,:), 'LineWidth', 1.8);
                         plot(ax, all_dl_nm{li,lj}, all_scat{li,lj},   '--', ...
                             'Color', colors(li,:), 'LineWidth', 1.2);

            g_fwhm = fwhm_g_nm(li, lj) * 1e3;
            s_fwhm = fwhm_s_nm(li, lj) * 1e3;
            if ~isnan(g_fwhm)
                leg_labels{li} = sprintf('L=%d \\mum  (G:%.0fpm  S:%.0fpm)', L_vals(li), g_fwhm, s_fwhm);
            else
                leg_labels{li} = sprintf('L=%d \\mum', L_vals(li));
            end
        end

        set(ax, 'YLim', [0 1.1]);  grid(ax, 'on');
        ylabel(ax, 'Norm. SHG Efficiency', 'FontSize', 10);
        title(ax, sprintf('\\alpha_0 = %d dB/cm  (SH = %d dB/cm)  —  %s  —  %s', ...
                           a0_sorted(lj), a0_sorted(lj)*ref_case.uv_loss_val, ...
                           ref_case.name, disp_note), 'FontSize', 8);
        legend(ax, h_leg, leg_labels, 'Location', 'northeast', 'FontSize', 8);
        text(0.02, 0.08, 'Solid = Guided,  Dashed = Scattered', ...
             'Units','normalized','FontSize',8,'Color',[0.35 0.35 0.35],'Parent',ax);

        if lj == n_loss_sel
            xlabel(ax, '\delta\lambda_{pump} (nm)  [detuning from perfect phase match;  \delta\lambda_{SHG} = \delta\lambda_{pump}/2]', 'FontSize', 10);
        end
    end

    fprintf('\nPhase Matching FWHM  (ref: %s)\n', ref_case.name);
    hdr = sprintf('%-10s', 'L (um)');
    for lj = 1:n_loss_sel
        hdr = [hdr sprintf('  | a0=%2ddB  G-FWHM(pm)  S-FWHM(pm)  S/G', a0_sorted(lj))]; %#ok
    end
    fprintf('%s\n%s\n', hdr, repmat('-', length(hdr)+2, 1));
    for li = 1:n_L
        row = sprintf('%-10d', L_vals(li));
        for lj = 1:n_loss_sel
            ratio = fwhm_s_nm(li,lj) / fwhm_g_nm(li,lj);
            row = [row sprintf('  |        %6.0f       %6.0f  %4.2f', ...
                               fwhm_g_nm(li,lj)*1e3, fwhm_s_nm(li,lj)*1e3, ratio)]; %#ok
        end
        fprintf('%s\n', row);
    end
    fprintf('%s\n', repmat('-', length(hdr)+2, 1));
end

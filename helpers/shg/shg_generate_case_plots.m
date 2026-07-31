function shg_generate_case_plots(case_struct, sweep_vals, guided, scat, blue_losses, param_label)
% Per-case figure: time-averaged power and normalized efficiency vs. sweep parameter.
    num_losses    = length(blue_losses);
    colors        = lines(num_losses);
    h_guided      = gobjects(num_losses, 1);
    legend_labels = cell(num_losses, 1);

    for c = 1:num_losses
        legend_labels{c} = sprintf('\\alpha_0=%d dB/cm', blue_losses(c));
    end

    fig_title = sprintf('Case: %s', case_struct.name);
    sub_text  = sprintf('P_{avg}=%.1fmW, DF=%.2g (P_{peak}=%.1fmW), UV loss: %s x%g', ...
                case_struct.Pp_avg_mW, case_struct.duty_factor, ...
                case_struct.Pp_avg_mW / case_struct.duty_factor, ...
                case_struct.uv_loss_mode, case_struct.uv_loss_val);

    figure('Color','w','Name',fig_title,'Position',[150, 150, 870, 680]);

    subplot(2,1,1); hold on;
    for c = 1:num_losses
        h_guided(c) = semilogy(sweep_vals, guided(:,c)*1e6, '-', 'Color', colors(c,:), 'LineWidth', 1.8);
                      semilogy(sweep_vals, scat(:,c)*1e6,   '--','Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on;
    ylabel('Time-Avg Power [\muW]');
    title({fig_title; sub_text}, 'FontSize', 9);
    legend(h_guided, legend_labels, 'Location', 'eastoutside', 'FontSize', 8);
    text(0.02, 0.92, 'Solid = Guided,  Dashed = Scattered', ...
         'Units','normalized','FontSize',8,'Color',[0.3 0.3 0.3]);

    subplot(2,1,2); hold on;
    Pp_avg_W = case_struct.Pp_avg_mW * 1e-3;
    h2 = gobjects(num_losses, 1);
    for c = 1:num_losses
        Eff_g = (guided(:,c) ./ Pp_avg_W^2) * 100;
        Eff_s = (scat(:,c)   ./ Pp_avg_W^2) * 100;
        h2(c) = semilogy(sweep_vals, Eff_g, '-', 'Color', colors(c,:), 'LineWidth', 1.8);
                semilogy(sweep_vals, Eff_s, '--','Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on;
    ylabel('Normalized Efficiency [%/W^2]'); xlabel(param_label);
    title('Time-Averaged Intrinsic Conversion Efficiency', 'FontSize', 9);
    legend(h2, legend_labels, 'Location', 'eastoutside', 'FontSize', 8);
end

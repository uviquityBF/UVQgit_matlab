function shg_generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses, target_L_um)
% SHG power vs. pump loss, at a fixed waveguide length (interpolated from the sweep).
    num_cases  = length(cases);
    num_losses = length(blue_losses);
    colors     = lines(num_cases);

    guided_at_L = zeros(num_cases, num_losses);
    scat_at_L   = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            guided_at_L(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L_um, 'linear', 'extrap');
            scat_at_L(k, c)   = interp1(sweep_vals, master_scat(:, c, k),   target_L_um, 'linear', 'extrap');
        end
    end

    [loss_sorted, sort_idx] = sort(blue_losses);
    guided_sorted = guided_at_L(:, sort_idx);
    scat_sorted   = scat_at_L(:,   sort_idx);

    fig_name = sprintf('SHG vs Loss at L = %d um', target_L_um);
    figure('Color','w','Name', fig_name,'Position',[200, 200, 820, 560]);
    hold on;

    h = gobjects(num_cases, 1);
    for k = 1:num_cases
        h(k) = semilogy(loss_sorted, guided_sorted(k,:)*1e6, '-o', ...
                        'Color', colors(k,:), 'LineWidth', 1.8, 'MarkerSize', 6, ...
                        'MarkerFaceColor', colors(k,:));
               semilogy(loss_sorted, scat_sorted(k,:)*1e6, '--o', ...
                        'Color', colors(k,:), 'LineWidth', 1.2, 'MarkerSize', 6, ...
                        'MarkerFaceColor', 'w');
    end

    set(gca, 'YScale', 'log'); grid on;
    xlabel('Pump Loss \alpha_0 (dB/cm)', 'FontSize', 11);
    ylabel('Time-Avg SHG Power (\muW)', 'FontSize', 11);
    title({sprintf('SHG Power vs Waveguide Loss  —  L = %d \\mum', target_L_um); ...
           'Solid filled = Guided,  Dashed open = Scattered'}, 'FontSize', 10);

    case_labels = cell(num_cases, 1);
    for k = 1:num_cases
        case_labels{k} = sprintf('%d: %s', k, cases(k).name);
    end
    legend(h, case_labels, 'Location', 'eastoutside', 'FontSize', 7);
end

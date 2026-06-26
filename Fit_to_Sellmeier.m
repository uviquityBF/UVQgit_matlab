function FitSellmeier_Range()
    %% 1. Setup and Load Data
    close all; clc;
    
    % --- CONFIGURATION ---
    defaultPath = 'C:\Users\brent\Downloads';
    if ~exist(defaultPath, 'dir'), startPath = pwd; else, startPath = defaultPath; end
    
    [file, path] = uigetfile(fullfile(startPath, '*.csv'), 'Select Dispersion CSV (WL, no, ne)');
    if isequal(file, 0), return; end
    
    dataRaw = importdata(fullfile(path, file));
    data = dataRaw.data;
    
    wl_raw = data(:,1); 
    no_raw = data(:,2); 
    ne_raw = data(:,3); 
    
    %% 2. Range Selection
    prompt = {'Start Wavelength (nm):', 'End Wavelength (nm):'};
    dlgtitle = 'Fitting Range Selection';
    definput = {'200', '600'};
    answer = inputdlg(prompt, dlgtitle, [1 40], definput);
    if isempty(answer), return; end
    
    wl_start = str2double(answer{1});
    wl_end = str2double(answer{2});
    
    mask = (wl_raw >= wl_start) & (wl_raw <= wl_end);
    wl = wl_raw(mask) / 1000; % Convert to um
    no = no_raw(mask);
    ne = ne_raw(mask);

    %% 3. Fit Logic
    sellmeier = @(x, wl) sqrt(1 + (x(1)*wl.^2)./(wl.^2 - x(2)) + ...
                                  (x(3)*wl.^2)./(wl.^2 - x(4)) + ...
                                  (x(5)*wl.^2)./(wl.^2 - x(6)));

    x0 = [1.0, 0.01, 0.8, 0.02, 1.0, 100]; 
    x0 = [1.0, 0.01, 0.8, 0.02, 1.0, 100]; 
    options = statset('MaxIter', 2000, 'TolFun', 1e-10);
    
    coeffs_o = nlinfit(wl, no, sellmeier, x0, options);
    coeffs_e = nlinfit(wl, ne, sellmeier, x0, options);

    %% 4. Visualization & Text Box
    wl_fit_um = linspace(min(wl), max(wl), 500);
    no_fit = sellmeier(coeffs_o, wl_fit_um);
    ne_fit = sellmeier(coeffs_e, wl_fit_um);

    hFig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 900 600]);
    hold on;
    
    plot(wl*1000, no, 'ko', 'MarkerSize', 4, 'DisplayName', 'n_o Data');
    plot(wl*1000, ne, 'ro', 'MarkerSize', 4, 'DisplayName', 'n_e Data');
    plot(wl_fit_um*1000, no_fit, 'k-', 'LineWidth', 2, 'DisplayName', 'n_o Fit');
    plot(wl_fit_um*1000, ne_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'n_e Fit');
    
    grid on; box on;
    xlabel('Wavelength (nm)'); ylabel('Refractive Index (n)');
    title(['Sellmeier Fit Results: ', file], 'Interpreter', 'none');
    legend('Location', 'northeast');

    % Construct the Results String for the Textbox
    % Note: B coefficients are dimensionless, C are in um^2
    resultsStr = {
        '\bfOrdinary Coefficients (n_o):', ...
        sprintf('  B1: %.5f, C1: %.5f', coeffs_o(1), coeffs_o(2)), ...
        sprintf('  B2: %.5f, C2: %.5f', coeffs_o(3), coeffs_o(4)), ...
        sprintf('  B3: %.5f, C3: %.5f', coeffs_o(5), coeffs_o(6)), ...
        '', ...
        '\bfExtraordinary Coefficients (n_e):', ...
        sprintf('  B1: %.5f, C1: %.5f', coeffs_e(1), coeffs_e(2)), ...
        sprintf('  B2: %.5f, C2: %.5f', coeffs_e(3), coeffs_e(4)), ...
        sprintf('  B3: %.5f, C3: %.5f', coeffs_e(5), coeffs_e(6)), ...
        '', ...
        sprintf('Range: %d - %d nm', wl_start, wl_end)
    };

    % Place the textbox in the bottom right (normalized coordinates)
    annotation('textbox', [0.55, 0.15, 0.33, 0.35], 'String', resultsStr, ...
        'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'Interpreter', 'tex', 'VerticalAlignment', 'middle');

    %% 5. Export to PNG
    [~, baseName, ~] = fileparts(file);
    outName = sprintf('%s_SellmeierFit_%d-%dnm.png', baseName, wl_start, wl_end);
    
    % Force the figure to render correctly before saving
    drawnow;
    saveas(hFig, fullfile(path, outName));
    
    fprintf('Fit Complete. Plot saved as: %s\n', outName);
end
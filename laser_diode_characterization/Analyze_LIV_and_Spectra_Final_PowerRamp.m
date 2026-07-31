function Analyze_LIV_and_Spectra_Final_PowerRamp()
    %% 1. Setup and Master File Selection
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\liv');
    close all; clear; clc;

    % --- DEFAULT CONFIGURATION ---
    powerScaleFactor = 1.437; 
    minThresholdFloormA = 20.0; % Search floor for "Flat" section (Current)
    rampStartPowermW = 15.0;    % Start of "Ramp" section (Power)
    defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\';
    
    if exist(defaultPath, 'dir'), startPath = defaultPath; else, startPath = pwd; end
    [masterFile, masterPath] = uigetfile(fullfile(startPath, '*.csv'), 'Select MASTER LIV CSV File');
    if isequal(masterFile, 0), return; end
    
    %% 2. Load Master LIV Data
    opts = detectImportOptions(fullfile(masterPath, masterFile));
    livTable = readtable(fullfile(masterPath, masterFile), opts);
    
    colNames = livTable.Properties.VariableNames;
    fileColIdx = width(livTable); 
    currIdx = find(contains(lower(colNames), 'actual_current'), 1);
    pwrIdx  = find(contains(lower(colNames), 'power_w'), 1);
    voltIdx = find(contains(lower(colNames), 'voltage'), 1);

    % Settings Prompt
    prompt = {'Start Wavelength (nm):', 'End Wavelength (nm):', 'Power Scaling:', ...
              'Search Floor for Flat Section (mA):', 'Ramp Start Power (mW):'};
    defAns = {'440', '460', num2str(powerScaleFactor), num2str(minThresholdFloormA), num2str(rampStartPowermW)};
    answer = inputdlg(prompt, 'Analysis Settings', [1 40], defAns);
    if isempty(answer), return; end
    
    wl_min = str2double(answer{1});
    wl_max = str2double(answer{2});
    powerScaleFactor = str2double(answer{3});
    minThresholdFloormA = str2double(answer{4});
    rampStartPowermW = str2double(answer{5});
    
    %% 3. Process Spectra
    numPoints = height(livTable);
    center_wls = NaN(numPoints, 1);
    hFig1 = figure('Name', 'Normalized Spectra Overlay', 'Color', 'w', 'Position', [100, 100, 1000, 600]); 
    hold on; grid on;

    for i = 1:numPoints
        rawSpecPath = livTable{i, fileColIdx}; 
        if iscell(rawSpecPath), rawSpecPath = rawSpecPath{1}; end
        if isempty(rawSpecPath) || (isnumeric(rawSpecPath) && isnan(rawSpecPath)), continue; end
        
        [~, name, ext] = fileparts(rawSpecPath);
        fullSpecPath = fullfile(masterPath, [name, ext]);
        
        if exist(fullSpecPath, 'file')
            dataRaw = importdata(fullSpecPath);
            if isstruct(dataRaw), rawMat = dataRaw.data; else, rawMat = dataRaw; end
            wl_raw = rawMat(:,1); int_raw = rawMat(:,2);
            
            [center_wls(i), maxVal, avg_noise, wl_seg, int_seg, is_valid] = spectral_centroid(wl_raw, int_raw, wl_min, wl_max);

            if is_valid
                int_norm = (int_seg - avg_noise) / (maxVal - avg_noise);
                hPlot = plot(wl_seg, int_norm, 'DisplayName', sprintf('%.1f mA', livTable{i, currIdx}));
                if numPoints > 15 && mod(i, 5) ~= 0 && i ~= numPoints && i ~= 1
                    set(get(get(hPlot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
            end
        end
    end
    title(['Spectral Overlay: ', masterFile], 'Interpreter', 'none');

    %% 4. Data Calculations & Dual-Linear Intersection (Power-Based)
    % --- Step A: Consolidate Duplicates ---
    [act_I, ~, idx_map] = unique(livTable{:, currIdx});
    act_P = accumarray(idx_map, livTable{:, pwrIdx} * 1000 * powerScaleFactor, [], @mean);
    act_V = accumarray(idx_map, livTable{:, voltIdx}, [], @mean);
    cwl_proc = accumarray(idx_map, center_wls, [], @(x) mean(x, 'omitnan'));

    % --- Step B: Dual-Linear Fit Algorithm ---
    % 1. Fit the "Floor" (Data points strictly below current search floor)
    floorMask = (act_I < minThresholdFloormA);
    if sum(floorMask) >= 2
        p_floor = polyfit(act_I(floorMask), act_P(floorMask), 1);
    else
        p_floor = [0, 0]; 
    end

    % 2. Fit the "Ramp" (Data points strictly above user Ramp Power threshold)
    rampMask = (act_P >= rampStartPowermW);
    if sum(rampMask) >= 2
        p_ramp = polyfit(act_I(rampMask), act_P(rampMask), 1);
        slope_eff_val = p_ramp(1); % mW/mA
    else
        p_ramp = [NaN, NaN];
        slope_eff_val = NaN;
    end

    % 3. Compute Intersection (m1*x + b1 = m2*x + b2)
    if ~isnan(p_ramp(1))
        % Intersection x = (b2 - b1) / (m1 - m2)
        I_threshold = (p_ramp(2) - p_floor(2)) / (p_floor(1) - p_ramp(1));
    else
        I_threshold = NaN;
    end

    % --- Step C: Extract Max Power & Spectral CWL ---
    [maxP_mW, maxIdxP] = max(act_P);
    cwl_at_maxP = cwl_proc(maxIdxP);

    % Wavelength Tuning Fit
    v_wl = ~isnan(cwl_proc);
    if any(v_wl)
        p_fit_wl = polyfit(act_I(v_wl), cwl_proc(v_wl), 1);
        fit_wl_plot = polyval(p_fit_wl, act_I);
    else
        p_fit_wl = [NaN, NaN];
        fit_wl_plot = NaN(size(act_I));
    end

    %% 5. Combined Plotting
    hFigCombo = figure('Name', 'Final Combined Results', 'Color', 'w', 'Position', [150, 150, 750, 900]);
    xLimits = [min(act_I) max(act_I)]; 
    
    subplot(2,1,1);
    yyaxis left
    plot(act_I, act_P, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    ylabel('True Power (mW)');
    hold on;
    % Mark the Intersection Threshold
    if ~isnan(I_threshold) && I_threshold > xLimits(1) && I_threshold < xLimits(2)
        y_cross = polyval(p_ramp, I_threshold);
        plot(I_threshold, y_cross, 'gx', 'MarkerSize', 14, 'LineWidth', 3, 'HandleVisibility','off');
    end
    
    yyaxis right
    plot(act_I, act_V, 'r-s', 'LineWidth', 1, 'MarkerFaceColor', 'r');
    ylabel('Voltage (V)');
    grid on;
    title(sprintf('%s\nLIV Summary (Scale: %.3f)', masterFile, powerScaleFactor), 'Interpreter', 'none');
    xlim(xLimits);
    
    livStatsStr = {sprintf('Max Power: %.2f mW', maxP_mW), ...
                   sprintf('Threshold (Dual-Fit): %.2f mA', I_threshold), ...
                   sprintf('Slope Eff: %.3f mW/mA', slope_eff_val)};
    annotation(hFigCombo, 'textbox', [0.15, 0.78, 0.28, 0.1], 'String', livStatsStr, ...
               'BackgroundColor', 'w', 'FontWeight', 'bold', 'FitBoxToText', 'on');
    
    subplot(2,1,2); hold on;
    if any(v_wl)
        plot(act_I(v_wl), cwl_proc(v_wl), 'ko', 'MarkerFaceColor', 'c', 'DisplayName', 'Measured');
        plot(act_I, fit_wl_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
        grid on; ylabel('Center WL (nm)'); xlabel('Actual Current (mA)');
        xlim(xLimits); legend('Location', 'best');
        
        specStatsStr = {sprintf('CWL @ Max P: %.3f nm', cwl_at_maxP), ...
                        sprintf('Tuning Rate: %.4e nm/mA', p_fit_wl(1))};
        annotation(hFigCombo, 'textbox', [0.40, 0.40, 0.20, 0.08], 'String', specStatsStr, ...
                   'BackgroundColor', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FitBoxToText', 'on');
    end

    %% 6. Final Export
    ts = datestr(now, 'yyyymmdd_HHMM');
    outXls = fullfile(masterPath, [strrep(masterFile, '.csv', ''), '_Final_', ts, '.xlsx']);
    cleanTable = table(act_I, act_V, act_P, cwl_proc, fit_wl_plot, ...
        'VariableNames', {'Current_mA', 'Voltage_V', 'True_Power_mW', 'Measured_WL_nm', 'Fit_WL_nm'});
    writetable(cleanTable, outXls);
    
    set(hFig1, 'PaperPositionMode', 'auto'); set(hFigCombo, 'PaperPositionMode', 'auto');
    print(hFig1, strrep(outXls, '.xlsx', '_Spectra.png'), '-dpng', '-r300');
    print(hFigCombo, strrep(outXls, '.xlsx', '_LIV.png'), '-dpng', '-r300');
    
    fprintf('\nExported: %s\nThreshold: %.2f mA\n', outXls, I_threshold);
end

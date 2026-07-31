function Analyze_LIV_Master()
    %% 1. Setup and File Selection
    close all; clear; clc;

    defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\nichia\NDBA116T (GainChips)\NichiaGC1_40C_auto_LIV_260407_spectra';
    if exist(defaultPath, 'dir'), startPath = defaultPath; else, startPath = pwd; end
    
    [file, path] = uigetfile(fullfile(startPath, '*.csv'), 'Select MASTER LIV CSV File');
    if isequal(file, 0), return; end
    
%% 2. Load Data (R2018b Compatible)
    opts = detectImportOptions(fullfile(path, file));
    
    % In R2018b, we skip 'VariableNamingRule'. 
    % MATLAB will automatically make headers "legal" variable names.
    livData = readtable(fullfile(path, file), opts);

 
    %% 3. Extract Columns
    % mapping based on your description
    I_set  = livData.set_current_mA;
    I_act  = livData.actual_current_mA;
    V_act  = livData.laser_actual_voltage_V;
    P_watt = livData.power_W;
    P_mw   = P_watt * 1000; % Convert to mW for standard LIV plotting

    %% 4. Visualization: Dual-Axis LIV Curve
    figure('Name', 'GainChip Master LIV', 'Color', 'w', 'Position', [100, 100, 900, 600]);
    
    % Light-Current (L-I)
    yyaxis left
    plot(I_act, P_mw, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    ylabel('Optical Power (mW)');
    set(gca, 'YColor', 'b');
    grid on;
    
    % Voltage-Current (V-I)
    yyaxis right
    plot(I_act, V_act, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    ylabel('Forward Voltage (V)');
    set(gca, 'YColor', 'r');
    
    xlabel('Actual Drive Current (mA)');
    title(['LIV Characteristic: ', file], 'Interpreter', 'none');

    %% 5. Calculate Slope Efficiency & Threshold
    % Slope Efficiency (mW/mA)
    % We use a moving window or simple diff
    dI = diff(I_act);
    dP = diff(P_mw);
    slopeEff = dP ./ dI;
    
    % Simple Threshold detection (Point of max change in slope)
    [~, maxD2Idx] = max(diff(slopeEff));
    threshold_I = I_act(maxD2Idx + 1);

    % Annotation Box
    statsStr = {sprintf('Threshold: %.2f mA', threshold_I), ...
                sprintf('Max Power: %.2f mW', max(P_mw)), ...
                sprintf('Max Voltage: %.2f V', max(V_act))};
            
    annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', statsStr, ...
               'BackgroundColor', 'w', 'FontWeight', 'bold', 'FitBoxToText', 'on');

    %% 6. Export Summary Table to Command Window
    % Selecting key columns for the summary
    Summary = table(I_set, I_act, V_act, P_mw, 'VariableNames', ...
        {'Set_mA', 'Actual_mA', 'Voltage_V', 'Power_mW'});
    disp(Summary);
end
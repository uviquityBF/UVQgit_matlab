%% Spiral Overhead Loss Analysis - R2018b Optimized
    % Specific Analysis for: (Run, Wafer, Subdie) combinations

function Historical_Loss_Data_Survey_analysis()
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_03_02_LabData\loss_recompile\'
    
    filename = '260306_SpiralOverhead_InitialSummaryStats - 260306_SpiralOverhead_InitialSummaryStats.csv';

    % 1. Load the Data
    data = readtable([dirs.d ,filename]);

    % 2. Clean Data (R^2 Filter)
    % Adjust this threshold as needed
    r2_threshold = 0.8;
    cleanData = data(data.R_2 > r2_threshold, :);

    % 3. Create the Composite Grouping Variable (Run + Wafer + Subdie)
    % In R2018b, we use string() for easy concatenation
    combos = string(cleanData.RunID) + " | " + ...
             string(cleanData.Wafer) + " | " + ...
             string(cleanData.SubdieNumber);

    cleanData.ComboGroup = categorical(combos);

    % 4. Prepare X-axis for Scatter Plot (Device ID)
    % Convert to categorical so MATLAB can plot string IDs on a numeric-like axis
    cleanData.DeviceID_Cat = categorical(cleanData.DeviceID);

    % 5. Visualization
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

    % --- Plot 1: Box Plot (Grouped by Combination) ---
    subplot(2,1,1);
    boxplot(cleanData.Loss_db_cm_, cleanData.ComboGroup);
    ylabel('Loss (db/cm)');
    title('Loss Distribution by [Run | Wafer | Subdie]');
    grid on;
    set(gca, 'TickLabelInterpreter', 'none');
    xtickangle(45); % Tilt labels so they don't overlap

    % --- Plot 2: Scatter Plot (Loss vs Device ID, Colored by Combination) ---
    subplot(2,1,2);
    % gscatter(X, Y, Group) creates a scatter plot with automatic coloring
    gscatter(cleanData.DeviceID_Cat, cleanData.Loss_db_cm_, cleanData.ComboGroup);

    xlabel('Device ID');
    ylabel('Loss (db/cm)');
    title('Device Loss Colored by [Run | Wafer | Subdie]');
    grid on;
    legend('Location', 'eastoutside', 'Interpreter', 'none');
    set(gca, 'TickLabelInterpreter', 'none');

    % Adjust layout to fit the legend
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2) pos(3)*0.85 pos(4)]);

    fprintf('Analysis Complete.\n');
    fprintf('Total unique combinations found: %d\n', numel(categories(cleanData.ComboGroup)));
end
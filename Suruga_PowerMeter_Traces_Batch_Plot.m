%% MATLAB Batch Optical Loss Data Processor (R2018b Compatible)
clear; clc; close all;

% =========================================================================
% CONFIGURATION
% =========================================================================
% Set your default directory path here (use \ for Windows, / for Mac/Linux)
defaultPath = 'C:\Your\Default\Folder\Path\Here'; 

% Fallback to current working directory if the default path doesn't exist
if ~exist(defaultPath, 'dir')
    defaultPath = pwd; 
end
% =========================================================================

% 1. Select specific files using uigetfile (MultiSelect enabled)
filePattern = fullfile(defaultPath, '*_PMlog.csv');
[fileNames, dataFolder] = uigetfile(filePattern, 'Select Optical Loss Logs', 'MultiSelect', 'on');

% Exit gracefully if user cancels
if isequal(fileNames, 0) || isequal(dataFolder, 0)
    error('No files selected. Exiting script.');
end

% Convert to cell array if only a single file was selected
if ischar(fileNames)
    fileNames = {fileNames};
end

numFiles = length(fileNames);

% 2. Initialize Figures for Plotting
% Figure 1: Absolute Calendar Time
fig1 = figure('Name', 'Optical Loss - Absolute Time', 'NumberTitle', 'off');
ax1 = axes(fig1); hold(ax1, 'on'); grid(ax1, 'on');
xlabel(ax1, 'Calendar Time');
ylabel(ax1, 'Loss / Ratio Out/In (dB)');
title(ax1, 'Optical Loss vs. Absolute Time');

% Figure 2: Relative Time (Starts at 0)
fig2 = figure('Name', 'Optical Loss - Relative Time', 'NumberTitle', 'off');
ax2 = axes(fig2); hold(ax2, 'on'); grid(ax2, 'on');
xlabel(ax2, 'Elapsed Time (seconds)'); % Adjust to minutes/hours if preferred
ylabel(ax2, 'Loss / Ratio Out/In (dB)');
title(ax2, 'Optical Loss vs. Relative Time (t_0 = 0)');

% 3. Preallocate arrays for the summary table
varNames = {'FabRun', 'Wafer', 'Subdie', 'DeviceID', 'Description', 'Min_Loss_dB', 'Max_Loss_dB', 'Mean_Loss_dB', 'Std_Dev'};

fabRunAll      = cell(numFiles, 1);
waferAll       = cell(numFiles, 1);
subdieAll      = cell(numFiles, 1);
deviceIDAll    = cell(numFiles, 1);
descriptionAll = cell(numFiles, 1);
minLossAll     = zeros(numFiles, 1);
maxLossAll     = zeros(numFiles, 1);
meanLossAll    = zeros(numFiles, 1);
stdDevAll      = zeros(numFiles, 1);

% 4. Loop through selected files
for k = 1:numFiles
    baseFileName = fileNames{k};
    fullFileName = fullfile(dataFolder, baseFileName);
    
    fprintf('Processing file %d of %d: %s\n', k, numFiles, baseFileName);
    
    %--- Parse Metadata from Filename ---
    cleanName = strrep(baseFileName, '.csv', '');
    nameParts = strsplit(cleanName, '_');
    
    if length(nameParts) >= 7
        fabRunAll{k}      = nameParts{2};
        waferAll{k}       = nameParts{3};
        subdieAll{k}      = nameParts{4};
        deviceIDAll{k}    = nameParts{5};
        descriptionAll{k} = strjoin(nameParts(6:end-1), '_');
    else
        fabRunAll{k}      = 'N/A';
        waferAll{k}       = 'N/A';
        subdieAll{k}      = 'N/A';
        deviceIDAll{k}    = 'N/A';
        descriptionAll{k} = 'Unknown Format';
    end
    
    %--- Read and Parse Data ---
    data = readtable(fullFileName);
    
    if iscell(data{:, 1}) || isstring(data{:, 1})
        timestamps = datetime(data{:, 1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    else
        timestamps = data{:, 1}; 
    end
    
    loss = data{:, 4};
    
    %--- Calculate Relative Time ---
    % Subtract the first timestamp element and convert the duration to seconds
    relativeTimeSec = seconds(timestamps - timestamps(1));
    
    %--- Plotting ---
    legendLabel = sprintf('W:%s, Dev:%s', waferAll{k}, deviceIDAll{k});
    
    % Plot 1: Absolute Calendar Time
    plot(ax1, timestamps, loss, 'LineWidth', 1.5, 'DisplayName', legendLabel);
    
    % Plot 2: Relative Time
    plot(ax2, relativeTimeSec, loss, 'LineWidth', 1.5, 'DisplayName', legendLabel);
    
    %--- Calculate Statistics ---
    minLossAll(k)  = min(loss);
    maxLossAll(k)  = max(loss);
    meanLossAll(k) = mean(loss);
    stdDevAll(k)   = std(loss);
end

% 5. Construct the final summary table
statsReport = table(fabRunAll, waferAll, subdieAll, deviceIDAll, descriptionAll, ...
                    minLossAll, maxLossAll, meanLossAll, stdDevAll, ...
                    'VariableNames', varNames);

% 6. Finalize Plot Features (Legends)
legend(ax1, 'Interpreter', 'none', 'Location', 'best');
legend(ax2, 'Interpreter', 'none', 'Location', 'best');

% 7. Display the Complete Statistics & Metadata Report
fprintf('\n========================================= STATISTICAL & METADATA REPORT =========================================\n');
disp(statsReport);
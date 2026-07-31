function Analyze_LIV_Spectra()
    %% 1. Setup and Path Configuration
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\liv');
    close all; clear; clc;

    defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\nichia\NDBA116T (GainChips)\NichiaGC1_40C_auto_LIV_260407_spectra';
    if exist(defaultPath, 'dir'), startPath = defaultPath; else, startPath = pwd; end
    
    %% 2. Select Files and Set Range
    [fileList, path] = uigetfile(fullfile(startPath, '*.csv'), ...
        'Select Spectra Files (Multi-Select)', 'MultiSelect', 'on');
    
    if isequal(fileList, 0), return; end
    if ischar(fileList), fileList = {fileList}; end 
    fileList = sort(fileList);

    prompt = {'Start Wavelength (nm):', 'End Wavelength (nm):'};
    dlgtitle = 'Adaptive Centroid Settings';
    definput = {'440', '460'}; 
    answer = inputdlg(prompt, dlgtitle, [1 40], definput);
    
    if isempty(answer), return; end
    wl_min = str2double(answer{1});
    wl_max = str2double(answer{2});

    %% 3. Pre-allocate Results
    numFiles = length(fileList);
    currents = zeros(numFiles, 1);
    center_wls = NaN(numFiles, 1); 
    max_ints = zeros(numFiles, 1);
    
    hFig1 = figure('Name', 'Normalized Spectra Overlay', 'Color', 'w');
    hold on; grid on;

    %% 4. Loop over files
    for i = 1:numFiles
        fileName = fileList{i};
        fullPath = fullfile(path, fileName);
        
        % Extract Current (Iset_...mA)
        tokens = regexp(fileName, 'Iset_([\d\.]+)mA', 'tokens');
        if ~isempty(tokens), currents(i) = str2double(tokens{1}{1}); else, currents(i) = NaN; end
        
        % Load Data
        dataRaw = importdata(fullPath);
        if isstruct(dataRaw), rawMat = dataRaw.data; else, rawMat = dataRaw; end
        wl_raw = rawMat(:,1);
        int_raw = rawMat(:,2);
        
        %% --- ADAPTIVE THRESHOLD & CENTROID (shared helper) ---
        [cwl, maxVal, avg_noise, wl_seg, int_seg, is_valid] = spectral_centroid(wl_raw, int_raw, wl_min, wl_max);

        if isempty(wl_seg), continue; end

        center_wls(i) = cwl;
        max_ints(i) = maxVal;

        if is_valid
            % Plotting: Normalize (0 to 1) for overlay
            int_norm = (int_seg - avg_noise) / (maxVal - avg_noise);
            plot(wl_seg, int_norm, 'DisplayName', sprintf('%.1f mA', currents(i)));
        end
    end

    %% 5. Finalize Visuals & Linear Fit
    figure(hFig1);
    xlabel('Wavelength (nm)'); ylabel('Normalized Intensity (Baseline Subtracted)');
    title('Adaptive Normalized Overlay');
    legend('show', 'Location', 'eastoutside');

    figure('Name', 'Center Wavelength vs. Current', 'Color', 'w');
    valid = ~isnan(center_wls) & ~isnan(currents);
    if any(valid)
        [plotX, sortIdx] = sort(currents(valid));
        tempWls = center_wls(valid);
        plotY = tempWls(sortIdx);
        
        % Perform Linear Fit (First order polynomial)
        p = polyfit(plotX, plotY, 1); 
        y_fit = polyval(p, plotX);
        tuningRate = p(1); % Slope: nm/mA
        
        % Plot Data and Fit
        hold on;
        plot(plotX, plotY, 'bo', 'MarkerFaceColor', 'c', 'DisplayName', 'Data');
        plot(plotX, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fit (%.4f nm/mA)', tuningRate));
        
        grid on; xlabel('Drive Current (mA)'); ylabel('Center Wavelength (nm)');
        title(['Wavelength Shift Fit: ', strrep(fileList{1}, '_', ' ')], 'Interpreter', 'none');
        legend('Location', 'best');
        
        % Stats annotation
        totalShift = max(plotY) - min(plotY);
        statsStr = {sprintf('Total Shift: %.4f nm', totalShift), ...
                    sprintf('Tuning Rate: %.4e nm/mA', tuningRate)};
        text(0.05, 0.95, statsStr, 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 10, ...
             'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    
    % Print Summary Table to Command Window
    fprintf('\n--- Analysis Summary ---\n');
    SummaryTable = table(currents(valid), max_ints(valid), center_wls(valid), ...
        'VariableNames', {'Current_mA', 'Peak_Intensity', 'Center_WL_nm'});
    disp(SummaryTable);
end

%% --- Helper Logic for R2018b ---
function val = ifdir(pth, tV, fV), if exist(pth, 'dir'), val = tV; else, val = fV; end; end
function val = ifstruct(s, tV, fV), if isstruct(s), val = tV; else, val = fV; end; end
function val = ifempty(v, tV, fV), if isempty(v), val = tV; else, val = fV; end; end

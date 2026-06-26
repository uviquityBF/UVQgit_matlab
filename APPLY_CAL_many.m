function APPLY_CAL_many()
    %% 1. Setup
    close all; clear; clc;
%     defaultPath = ...
%     'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Tools\_SYSTEM-2\UV Response Calibration of Detection [cts uJ]\process Cal in matlab\Sys2_Calibration_2026_03_17 (50um slit)';
    defaultPath = ...
    'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Tools\_SYSTEM-2\UV Response Calibration of Detection [cts uJ]\process Cal in matlab\Sys2_Calibration_2026_03_19 (50um slit)';

%     defaultPath = ...
%     'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Tools\_SYSTEM-2\UV Response Calibration of Detection [cts uJ]\process Cal in matlab\Sys2_Calibratin_2026_05_05 (50um slit)'
    if ~exist(defaultPath, 'dir'), startPath = pwd; else, startPath = defaultPath; end
    
    %% 2. Read Inputs
    [file1, path1] = uigetfile(fullfile(startPath, '*.csv'), 'Select RESPONSE of Spectrometer ( CSV )');
    if isequal(file1, 0), return; end
    raw1 = importdata(fullfile(path1, file1));
    A_Resp = raw1.data; 
    
    [file2, path2] = uigetfile(fullfile(path1, '*.txt'), 'Select Fiber10 Spectrum (QEPro TXT)');
    if isequal(file2, 0), return; end
    fiberData = readQEPro_SpectrumFile(fullfile(path2, file2));
    
    [fileList, path3] = uigetfile(fullfile(path2, '*.txt'), 'Select ALL DUT Spectra', 'MultiSelect', 'on');
    if isequal(fileList, 0), return; end
    if ischar(fileList), fileList = {fileList}; end

    %% 3. Initialize Results Table and Figures
    firstData = readQEPro_SpectrumFile(fullfile(path3, fileList{1}));
    wl_master = firstData.A(:, 1);
    resultsTable = table(wl_master, 'VariableNames', {'Wavelength_nm'});
    
    
    hFig1 = figure('Name', 'Responsivity Comparison'); 
    ax1 = axes('Parent', hFig1, 'YScale', 'log'); hold(ax1, 'on');
    
    hFig2 = figure('Name', 'Optics Transmission'); 
    ax2 = axes('Parent', hFig2, 'YScale', 'log'); hold(ax2, 'on'); 

    %% 3b. Get User Input for Total Power
%     % Power constant and string for labeling
%     PowerTot_uW_Fiber10 = 0.011; 
%     pwrStr = strrep(num2str(PowerTot_uW_Fiber10), '.', 'p'); % e.g., "0p011"

    % Define the default value
    defaultPower = 0.011; 
    
    % Setup the dialog box
    prompt = {'Enter Total Spectral Power from Fiber10 [uW]:'};
    dlgtitle = 'Power Configuration';
    dims = [1 50];
    definput = {num2str(defaultPower)};
    
    % Show the dialog
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    % If user cancels, stop the script. Otherwise, convert to double.
    if isempty(answer)
        fprintf('User cancelled. Exiting...\n');
        return;
    else
        PowerTot_uW_Fiber10 = str2double(answer{1});
    end

    % Update the Power String for filenames (e.g., 0.011 -> 0p011)
    pwrStr = strrep(num2str(PowerTot_uW_Fiber10), '.', 'p');    
    
    
%% 4. Loop Through Each DUT File
    for i = 1:length(fileList)
        
        thisFile = fileList{i};
        
        % Extract Z value for Column Header
        z_match = regexp(thisFile, 'Z=([\d\.]+)', 'tokens');
        if ~isempty(z_match)
            colName = ['Z_', strrep(z_match{1}{1}, '.', 'p')];
            displayName = ['Z = ', z_match{1}{1}];
        else
            colName = sprintf('Sample_%02d', i);
            displayName = strrep(thisFile, '_', ' ');
        end
        
        % Load and Calibrate
        dutData = readQEPro_SpectrumFile(fullfile(path3, thisFile));
        this_counts = dutData.A(:, 2);
        this_Tint = dutData.Tint;

        % 2. Calculate Pixel Width (Jacobian)
        % This accounts for non-linear dispersion (nm per pixel)
        dl = zeros(size(wl_master));
        dl(1:end-1) = diff(wl_master); 
        dl(end) = dl(end-1); % Fill last element to maintain array size
        
        % 3. Interpolate Static Cal Data
        resp_interp    = interp1(A_Resp(:,1), A_Resp(:,2), wl_master, 'linear', 'extrap');
        fiber10_interp = interp1(fiberData.A(:,1), fiberData.A(:,2), wl_master, 'linear', 'extrap');

        % 4. Calibration Logic
        dut_cps = this_counts ./ this_Tint;
        
        % Masking to avoid noise at the edges (<202nm) and NaNs
        isValid = ~isnan(resp_interp) & resp_interp>0 & ~isnan(fiber10_interp) & (wl_master > 202);

        % 5. Distribute Total Power [uW] across the spectrum
        % We weight the spectrum by pixel width (nm) to get a "Power per nm" distribution
        CAL_weighted = zeros(size(wl_master));
%         CAL_weighted(isValid) = fiber10_interp(isValid) .* dl(isValid);
%        CAL_weighted(isValid) = fiber10_interp(isValid) .* dl(isValid) ./ (resp_interp(isValid)+eps);
        
        total_area = sum(CAL_weighted(isValid)); % Total "Area" in Count-nanometers
        
        % Normalize and scale by the measured Total Power 
        % This gives us [uW / sample] for every wavelength bin
        CAL_Spectrum_power = (CAL_weighted ./ total_area) * PowerTot_uW_Fiber10;
        CAL_Spectrum_power_inband = CAL_Spectrum_power .* ((295 <= wl_master & wl_master <=305)) * 0.4
        

        % 5. Calculate Final Responsivity [cts / uJ]
        Response_DUT = zeros(size(wl_master));
        canCalc = isValid & (CAL_Spectrum_power > 0);
        Response_DUT(canCalc) = dut_cps(canCalc) ./ CAL_Spectrum_power(canCalc);  %[cts / uJ]
        
        % Store and Plot...
        resultsTable.(colName) = Response_DUT;
        semilogy(ax1, wl_master, Response_DUT, 'LineWidth', 1.2, 'DisplayName', displayName);
        semilogy(ax2, wl_master, Response_DUT ./ (resp_interp + eps), 'LineWidth', 1.2, 'DisplayName', displayName);
    end
    
%     figure;  plot(wl_master, CAL_weighted);
    
    %% 5. Overlay Reference and Finalize Visuals
    % Overlay Spectrometer Responsivity on Figure 1
    semilogy(ax1, A_Resp(:,1), A_Resp(:,2), 'k--', 'LineWidth', 2.5, 'DisplayName', 'Spec Responsivity (Ref)');

    % Format Figure 1
    grid(ax1, 'on'); grid(ax1, 'minor');
    ylabel(ax1, 'Response [cts/uJ]'); xlabel(ax1, 'Wavelength (nm)');
    xlim(ax1, [200, 400]); ylim(ax1, [1e6, 1e9]);
    title(ax1, {['Spectral Responsivity Comparison (P=', num2str(PowerTot_uW_Fiber10), ' uW)']; 'Reference vs DUTs'});
    legend(ax1, 'show', 'Location', 'eastoutside');

    % Format Figure 2
    grid(ax2, 'on'); grid(ax2, 'minor');
    ylabel(ax2, 'Transmission Ratio (DUT / Ref)'); xlabel(ax2, 'Wavelength (nm)');
    xlim(ax2, [200, 400]); ylim(ax2, [1e-3, 2]); % Adjusted for log scale visibility
    title(ax2, ['Relative Transmission of Optics (P=', num2str(PowerTot_uW_Fiber10), ' uW)']);
    legend(ax2, 'show', 'Location', 'eastoutside');
    %% 5.5 TEMP OUTPUTS
    wl2check = 223;
    [dummy,isel] = min(abs(wl_master-wl2check));
    disp(['CAL spectrum fraction at ',num2str(wl2check),':   ',num2str( CAL_Spectrum_power(isel)/PowerTot_uW_Fiber10 ),'[1/bin]']);
    
    tmp = resultsTable{isel,:};
    DUT_cps_at_wl2check = max(tmp)*CAL_Spectrum_power(isel);
    disp(['DUT_cps_at_wl2check: ',num2str(DUT_cps_at_wl2check),'cps'])
   
    pause(1);
    
    %% 6. Export
    outPath = path1; % Saving to the DUT folder
    csvName = sprintf('MultiBatch_Results_P%s_uW.csv', pwrStr);
    fig1Name = sprintf('Plot_Responsivity_P%s_uW.png', pwrStr);
    fig2Name = sprintf('Plot_Transmission_P%s_uW.png', pwrStr);
    
    writetable(resultsTable, fullfile(outPath, csvName));
    saveas(hFig1, fullfile(outPath, fig1Name));
    saveas(hFig2, fullfile(outPath, fig2Name));
    
    fprintf('Processed %d files. CSV and PNGs saved with Power Tag: %s\n', length(fileList), pwrStr);
end

%% --- Helper Functions ---
function dataObj = readQEPro_SpectrumFile(fullPath)
    dataObj.A = []; dataObj.Tint = NaN;
    fid = fopen(fullPath, 'r'); if fid == -1, return; end
    lineIdx = 0; dataStartLine = 0;
    while ~feof(fid) && lineIdx < 100
        line = fgetl(fid); lineIdx = lineIdx + 1;
        if contains(line, '(sec):')
            parts = strsplit(line, ':'); 
            dataObj.Tint = str2double(parts{end}); 
        end
        if ~isempty(line) && (isstrprop(line(1), 'digit') || line(1) == '-')
            dataStartLine = lineIdx; break; 
        end
    end
    fclose(fid);
    if dataStartLine > 0
        try
            dataObj.A = dlmread(fullPath, '\t', dataStartLine - 1, 0);
        catch
            dataObj.A = manual_parse(fullPath, dataStartLine);
        end
    end
end

function A = manual_parse(fullPath, startLine)
    fid = fopen(fullPath, 'r'); 
    for i = 1:startLine-1, fgetl(fid); end 
    raw = textscan(fid, '%f %f'); 
    fclose(fid); A = [raw{1}, raw{2}];
end

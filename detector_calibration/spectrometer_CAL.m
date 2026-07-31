function MAIN()
    %% 1. Setup and Path Configuration
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\spectra');
    close all; clear; clc;

%     defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Tools\UV Collection Calibration [cts uJ]\process Cal in matlab';
%     defaultPath  = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Tools\UV Collection Calibration [cts uJ]\process Cal in matlab\Maya_Calibration_2026_03_26';
%     defaultPath  = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Components\Spectrometer -- MAYA\UV Responsivity Calibration\process Cal in Matlab\Maya_Calibration_2026_03_11 (50um)';
%     defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Components\Spectrometer -- MAYA\UV Responsivity Calibration\process Cal in Matlab\Maya_Calibration_2026_03_26 (10um, 200um, metal)\';
      defaultPath  = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Components\Spectrometer -- QE_Pro\UV Responsivity Calibration\process Cal in Matlab\QEPro_Calibration_2026_03_11 (50um)';
%       defaultPath  = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Components\Spectrometer -- QE_Pro\UV Responsivity Calibration\process Cal in Matlab\QEPro_Calibration_2026_03_24 (50um)'
%     defaultPath = 'G:\Shared drives\Corp Main\Engineering\LAB Equipment\Components\Spectrometer -- QE_Pro\UV Responsivity Calibration\process Cal in Matlab\QEPro_Calibration_2026_05_04 (50um)';
    if ~exist(defaultPath, 'dir')
        startPath = pwd; 
    else
        startPath = defaultPath;
    end
    
    %% 2. Read Inputs
    currentPath = startPath; 

    [file1, path1] = uigetfile(fullfile(currentPath, '*.csv'), 'Select SL3-CAL Reference Spectrum');
    if isequal(file1, 0), error('no reference spectrum provided. Ending.'); return; end
    raw1 = importdata(fullfile(path1, file1));
    A1 = raw1.data; 
    
    [file2, path2] = uigetfile(fullfile(path1, '*.csv'), 'Select CR2 Transmission');
    if isequal(file2, 0)
        A2 = [[100:1:1000]',ones(901,1)];
        disp('No CR2 Transmission selected.  Setting T=100%'); 
    else
        raw2 = importdata(fullfile(path2, file2));
        A2 = raw2.data;
    end

    [file3, path3] = uigetfile(fullfile(path2, '*.csv'), 'Select Fiber C Transmission');
    if isequal(file3, 0)
        A3 = [[100:1:1000]',ones(901,1)];
       disp('No FiberC Transmission selected.  Setting T=100%');
    else
        raw3 = importdata(fullfile(path3, file3));
        A3 = raw3.data;
    end
    
    % Load DUT (using your specific QEPro helper)
    [file4, path4] = uigetfile(fullfile(path3, '*.txt'), 'Select DUT Spectrum (QEPro)');
    if isequal(file4, 0), return; end
    dutData = read_spectrum_file(fullfile(path4, file4), 'oceanoptics');
    IntegrationTime = dutData.Tint;

    %% 3. Process Input Data
    wl_master = dutData.wl;
    dut_counts = dutData.counts;
    
    ref_spectrum_interp = interp1(A1(:,1), A1(:,2), wl_master, 'linear', 'extrap');
    cr2_trans_interp    = interp1(A2(:,1), A2(:,2), wl_master, 'linear', 'extrap');
    fiber_trans_interp  = interp1(A3(:,1), A3(:,2), wl_master, 'linear', 'extrap');

    %% 3b. Visualize Input Data (Figure 1)
    hFig1 = figure('Name', 'Input Data Check', 'Visible', 'on');
    subplot(3,1,1);
    plot(wl_master, ref_spectrum_interp, 'k');
    ylabel('[W/m2]'); title('INPUT: Reference Spectra'); grid on;

    subplot(3,1,2);
    plot(wl_master, cr2_trans_interp, 'r'); hold on;
    plot(wl_master, fiber_trans_interp, 'b');
    ylabel('Transmission'); title('INPUT: Component Transmission'); grid on; ylim([0,1.1]);
    
    subplot(3,1,3);
    plot(wl_master, dut_counts, 'g');
    ylabel('Counts'); title('INPUT: DUT Raw Spectrum'); xlabel('Wavelength (nm)'); grid on;
 
    %% 3b --- Total Spectral Power Measured
	% Define the default value
    defaultPower = 0.066; 
    
    % Setup the dialog box
    prompt = {'Enter Total Spectral Power from FiberC [uW]:'};
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
        PowerTot_uW_SL3CAL_CR2_FiberC = str2double(answer{1});
    end

    %% 4. Calibration Math
    dut_cps = dut_counts ./ IntegrationTime;

    isValid = ~isnan(ref_spectrum_interp) & ...
              ~isnan(cr2_trans_interp) & ...
              ~isnan(fiber_trans_interp) & ...
              (wl_master > 202);

    Response_DUT_cts_per_uJ = zeros(size(wl_master));
    CAL_Spectrum_raw = zeros(size(wl_master));

    CAL_Spectrum_raw(isValid) = ref_spectrum_interp(isValid) .* ...
                                 cr2_trans_interp(isValid) .* ...
                                 fiber_trans_interp(isValid);
                             
    total_raw_sum = sum(CAL_Spectrum_raw(isValid));
    CAL_Spectrum_normalized = CAL_Spectrum_raw ./ total_raw_sum;
    CAL_Spectrum_power = CAL_Spectrum_normalized * PowerTot_uW_SL3CAL_CR2_FiberC;

    canCalculateR = isValid & (CAL_Spectrum_power > 0);
    Response_DUT_cts_per_uJ(canCalculateR) = dut_cps(canCalculateR) ./ CAL_Spectrum_power(canCalculateR);

    %% 5. Visualization (Figure 2)
    hFig2 = figure('Name', 'Final Calibration Result');
    semilogy(wl_master, Response_DUT_cts_per_uJ, 'LineWidth', 1.5); 
    grid on;
    ylabel('DUT Response [cts/uJ]'); xlabel('Wavelength (nm)');
    ylim([1e7, 1e10]); xlim([200, 700]);
    title({['Responsivity: ', file4];['assumed ref power [uW]:',num2str(PowerTot_uW_SL3CAL_CR2_FiberC),'uW']}, 'Interpreter', 'none');

    %% 5.5 TEMP OUTPUTS
    wl2check = 228;
    [dummy,isel] = min(abs(wl_master-wl2check));
    disp(['CAL_Spectrum_normalized: ',num2str(CAL_Spectrum_normalized(isel)),' [a.u. unitless]'])
    disp(['QE Pro measured signal: ',num2str(dut_cps(isel)),' [cps]'])
     
    
    pause(1);
    
    %% 6. Automatic Export
    % Use the DUT filename (file4) to create output names
    [~, baseName, ~] = fileparts(file4);
    
    % Save CSV
    outputTable = table(wl_master, Response_DUT_cts_per_uJ, ...
        'VariableNames', {'Wavelength_nm', 'Response_cts_per_uJ'});
    csvOutName = fullfile(path4, [baseName, '_Calibration.csv']);
    writetable(outputTable, csvOutName);
    
    % Save PNGs
    saveas(hFig1, fullfile(path4, [baseName, '_InputCheck.png']));
    saveas(hFig2, fullfile(path4, [baseName, '_FinalResponse.png']));
    
    fprintf('Analysis Complete.\nFiles saved to: %s\n', path4);
end
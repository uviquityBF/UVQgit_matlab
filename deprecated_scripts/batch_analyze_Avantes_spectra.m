%% Open Avantes Spectra (ASCII) with Gaussian and Centroid Analysis

function batch_analyze_spectra_centroid()
    close all
    % --- Setup Paths ---
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM135\SD48\GD2\260408_50um_TiSaph_WindowSweep\Avavantes_PumpSpectrum';
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\nichia\NDBA116T (GainChips)\GainSpectra';
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_04_06_LabData\ECL with nichia GC2\'
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136\SD49\GD3\251124\SHG Search\id38.3_50um_Maya';
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_04_27_LabData\DLPro Tuning';
    default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_04_27_LabData\Avantes -v- OSA -v- Wavemeter\Avantes';
    type = 'Avantes'; %'OceanOptics'; %'Avantes'; %'OceanOptics'

    if ~exist(default_dir, 'dir'); default_dir = pwd; end

    % 1. Multi-select files
    [filenames, start_path] = uigetfile({'*.txt', 'Avantes Spectra (*.txt)'}, ...
        'Select Spectra Files', default_dir, 'MultiSelect', 'on');
    
    if isequal(filenames, 0); return; end
    if ischar(filenames); filenames = {filenames}; end 
    
    num_files = length(filenames);
    target_int_time = 100; 
    
    % Initialize containers
    all_WL = [];
    all_normalized = [];
    all_time_scaled = [];
    labels = cell(1, num_files);
    fit_center_wl = zeros(1, num_files);
    centroid_wl = zeros(1, num_files); % New container for Centroid

    fprintf('\n--- Processing Results (Pump Spectrum Analysis) ---\n');
    fprintf('%-40s | %-12s | %-12s\n', 'Filename', 'Gaussian nm', 'Centroid nm');
    fprintf('%s\n', repmat('-', 1, 75));

    for k = 1:num_files
        full_path = fullfile(start_path, filenames{k});
        labels{k} = strrep(filenames{k}, '_', ' ');

        % --- Extract Integration Time ---
        int_time = 1; 
        fid = fopen(full_path, 'r');
        for i = 1:9
            line_str = fgetl(fid);
            if ~isempty(strfind(lower(line_str), 'integration time'))
                nums = sscanf(line_str, '%*[^0-9.]%f'); 
                if ~isempty(nums); int_time = nums(1); end
            end
        end
        fclose(fid);

        % --- Read Numeric Data ---
        switch type
            case 'Avantes'
                data = dlmread(full_path, ';', 9,0);
            case 'OceanOptics'
                raw = importdata(full_path, '\t', 17);
                data  = raw.data;
        end
        
        WL = data(:,1);
        Intensity = data(:,2);
        if k == 1; all_WL = WL; end 

        % 2. Gaussian Fit (430nm to 470nm)
        idx = find(WL >= 430 & WL <= 470);
        if ~isempty(idx) && length(idx) > 3
            x = WL(idx);
            y = Intensity(idx);
            y_fit = y;
            y_fit(y_fit <= 0) = min(y_fit(y_fit > 0)); 
            logY = log(y_fit);
            coeffs = polyfit(x, logY, 2);
            fit_center_wl(k) = -coeffs(2) / (2 * coeffs(1));
        else
            fit_center_wl(k) = NaN;
        end

        % --- ***CENTROID CALCULATION HERE*** ---
        % A. Define Baseline (200-300nm) to characterize noise floor
        baseMask = (WL >= 200) & (WL <= 300);
        if any(baseMask)
            avg_noise = mean(Intensity(baseMask));
            std_noise = std(Intensity(baseMask));
            auto_thresh = avg_noise + (10 * std_noise);
        else
            avg_noise = 0; auto_thresh = 100; % Fallback
        end

        % B. Use 50% Threshold logic within the peak segment
        if ~isempty(idx)
            wl_seg = WL(idx);
            int_seg = Intensity(idx);
            [maxVal, ~] = max(int_seg);
            
            if maxVal > auto_thresh
                % Threshold relative to baseline
                cutOff = avg_noise + (maxVal - avg_noise) * 0.5;
                peakMask = int_seg >= cutOff;
                
                % Weighted average (Centroid)
                wl_fc = wl_seg(peakMask);
                int_fc = int_seg(peakMask);
                centroid_wl(k) = sum(wl_fc .* int_fc) / sum(int_fc);
            else
                centroid_wl(k) = NaN;
            end
        else
            centroid_wl(k) = NaN;
        end
        
        % Print both results for comparison
        fprintf('%-40s | %-12.2f | %-12.2f\n', filenames{k}, fit_center_wl(k), centroid_wl(k));

        % Rescaling logic
        max_val = max(Intensity);
        if max_val == 0; max_val = 1; end 
        all_normalized(:, k) = Intensity / max_val;
        all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
    end

    % --- Plotting ---
    plot_spectra(all_WL, all_normalized, 'Normalized (Max=1)', labels);
    plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);

    % 3. Optional CSV Output
    handle_csv_export(all_WL, all_normalized, all_time_scaled, filenames);
end

% --- Sub-function: CSV Handling ---
function handle_csv_export(wl, norm_data, time_data, filenames)
    answer = questdlg('Export processed data to two CSV files?', 'Export Data', 'Yes', 'No', 'No');
    if strcmp(answer, 'Yes')
        [save_name, save_path] = uiputfile('Spectra_MaxUnity.csv', 'Save Normalized Data As');
        if ~isequal(save_name, 0)
            export_csv_with_headers(fullfile(save_path, save_name), wl, norm_data, filenames);
            time_scaled_name = strrep(save_name, 'MaxUnity', 'TimeCorrected');
            if strcmp(time_scaled_name, save_name); time_scaled_name = ['TimeCorrected_', save_name]; end
            export_csv_with_headers(fullfile(save_path, time_scaled_name), wl, time_data, filenames);
        end
    end
end

function export_csv_with_headers(full_file_path, wl, data_matrix, headers)
    fid = fopen(full_file_path, 'w');
    fprintf(fid, 'Wavelength_nm');
    for i = 1:length(headers); fprintf(fid, ',%s', headers{i}); end
    fprintf(fid, '\n');
    fclose(fid);
    dlmwrite(full_file_path, [wl, data_matrix], '-append', 'delimiter', ',');
end

function plot_spectra(x, y, title_str, legs)
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
    plot(x, y, 'LineWidth', 2); 
    grid on; xlim([435,455]);
    xlabel('Wavelength (nm)'); ylabel('Intensity');
    title(title_str);
    legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
end










% function batch_analyze_spectra_gaussian()
%     close all
%     % --- Setup Paths ---
% %     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_23_LabData\bluglass_lot2_GC1_vbg_improved_collimation';
%     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM135\SD48\GD2\260408_50um_TiSaph_WindowSweep\Avavantes_PumpSpectrum';
%     if ~exist(default_dir, 'dir'); default_dir = pwd; end
% 
%     % 1. Multi-select files
%     [filenames, start_path] = uigetfile({'*.txt', 'Avantes Spectra (*.txt)'}, ...
%         'Select Spectra Files', default_dir, 'MultiSelect', 'on');
%     
%     if isequal(filenames, 0); return; end
%     if ischar(filenames); filenames = {filenames}; end 
%     
%     num_files = length(filenames);
%     target_int_time = 100; 
%     
%     % Initialize containers
%     all_WL = [];
%     all_normalized = [];
%     all_time_scaled = [];
%     labels = cell(1, num_files);
%     fit_center_wl = zeros(1, num_files);
% 
%     fprintf('\n--- Processing Results (Gaussian Fit 430-470nm) ---\n');
%     fprintf('%-40s | %-20s\n', 'Filename', 'Peak WL (Gaussian)');
%     fprintf('%s\n', repmat('-', 1, 65));
% 
%     for k = 1:num_files
%         full_path = fullfile(start_path, filenames{k});
%         labels{k} = strrep(filenames{k}, '_', ' ');
% 
%         % --- Extract Integration Time ---
%         int_time = 1; 
%         fid = fopen(full_path, 'r');
%         for i = 1:9
%             line_str = fgetl(fid);
%             if ~isempty(strfind(lower(line_str), 'integration time'))
%                 nums = sscanf(line_str, '%*[^0-9.]%f'); 
%                 if ~isempty(nums); int_time = nums(1); end
%             end
%         end
%         fclose(fid);
% 
%         % --- Read Numeric Data ---
%         data = dlmread(full_path, ';', 9, 0);
%         WL = data(:,1);
%         Intensity = data(:,2);
%         if k == 1; all_WL = WL; end 
% 
%         % 2. Gaussian Fit (430nm to 470nm)
%         idx = find(WL >= 430 & WL <= 470);
%         if ~isempty(idx) && length(idx) > 3
%             x = WL(idx);
%             y = Intensity(idx);
%             
%             % Ensure no zeros/negatives for the log transform
%             y(y <= 0) = min(y(y > 0)); 
%             logY = log(y);
%             
%             % Fit a 2nd order polynomial: log(y) = ax^2 + bx + c
%             coeffs = polyfit(x, logY, 2);
%             a = coeffs(1);
%             b = coeffs(2);
%             
%             % The peak of the parabola is at -b / 2a
%             fit_center_wl(k) = -b / (2 * a);
%         else
%             fit_center_wl(k) = NaN;
%         end
% 
%         %***CENTROID CALCULATION HERE***
%         
%         fprintf('%-40s | %-6.2f nm\n', filenames{k}, fit_center_wl(k));
% 
%         % Rescaling logic
%         max_val = max(Intensity);
%         if max_val == 0; max_val = 1; end 
%         all_normalized(:, k) = Intensity / max_val;
%         all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
%     end
% 
%     % --- Plotting ---
%     plot_spectra(all_WL, all_normalized, 'Normalized (Max=1)', labels);
%     plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);
% 
%     % 3. Optional CSV Output
%     handle_csv_export(all_WL, all_normalized, all_time_scaled, filenames);
% end
% 
% % --- Sub-function: CSV Handling ---
% function handle_csv_export(wl, norm_data, time_data, filenames)
%     answer = questdlg('Export processed data to two CSV files?', 'Export Data', 'Yes', 'No', 'No');
%     if strcmp(answer, 'Yes')
%         [save_name, save_path] = uiputfile('Spectra_MaxUnity.csv', 'Save Normalized Data As');
%         if ~isequal(save_name, 0)
%             export_csv_with_headers(fullfile(save_path, save_name), wl, norm_data, filenames);
%             
%             time_scaled_name = strrep(save_name, 'MaxUnity', 'TimeCorrected');
%             if strcmp(time_scaled_name, save_name); time_scaled_name = ['TimeCorrected_', save_name]; end
%             
% %             [save_name2, save_path2] = uiputfile(time_scaled_name, 'Save Time-Corrected Data As');
%             save_name2 = time_scaled_name;
%             if ~isequal(save_name2, 0)
%                 export_csv_with_headers(fullfile(save_path, save_name2), wl, time_data, filenames);
%             end
%         end
%     end
% end
% 
% function export_csv_with_headers(full_file_path, wl, data_matrix, headers)
%     fid = fopen(full_file_path, 'w');
%     fprintf(fid, 'Wavelength_nm');
%     for i = 1:length(headers); fprintf(fid, ',%s', headers{i}); end
%     fprintf(fid, '\n');
%     fclose(fid);
%     dlmwrite(full_file_path, [wl, data_matrix], '-append', 'delimiter', ',');
% end
% 
% function plot_spectra(x, y, title_str, legs)
%     figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
%     plot(x, y, 'LineWidth', 2); 
%     grid on; xlim([435,455]);
%     xlabel('Wavelength (nm)'); ylabel('Intensity');
%     title(title_str);
%     legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
% end
% 
% 
% 










% function batch_analyze_spectra_final()
%     % --- Setup Paths ---
%     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_16_LabData';
%     if ~exist(default_dir, 'dir'); default_dir = pwd; end
% 
%     % 1. Multi-select files
%     [filenames, start_path] = uigetfile({'*.txt', 'Avantes Spectra (*.txt)'}, ...
%         'Select Spectra Files', default_dir, 'MultiSelect', 'on');
%     
%     if isequal(filenames, 0); return; end
%     if ischar(filenames); filenames = {filenames}; end 
%     
%     num_files = length(filenames);
%     target_int_time = 100; 
%     
%     % Initialize containers
%     all_WL = [];
%     all_normalized = [];
%     all_time_scaled = [];
%     labels = cell(1, num_files);
%     center_wavelengths = zeros(1, num_files);
% 
%     fprintf('\n--- Processing Results ---\n');
%     fprintf('%-40s | %-20s\n', 'Filename', 'Center WL (420-470nm)');
%     fprintf('%s\n', repmat('-', 1, 65));
% 
%     for k = 1:num_files
%         full_path = fullfile(start_path, filenames{k});
%         labels{k} = strrep(filenames{k}, '_', ' ');
% 
%         % --- Extract Integration Time ---
%         int_time = 1; 
%         fid = fopen(full_path, 'r');
%         for i = 1:9
%             line_str = fgetl(fid);
%             if ~isempty(strfind(lower(line_str), 'integration time'))
%                 nums = sscanf(line_str, '%*[^0-9.]%f'); 
%                 if ~isempty(nums); int_time = nums(1); end
%             end
%         end
%         fclose(fid);
% 
%         % --- Read Numeric Data ---
%         data = dlmread(full_path, ';', 9, 0);
%         WL = data(:,1);
%         Intensity = data(:,2);
%         if k == 1; all_WL = WL; end 
% 
%         % 2. Center of Mass (420nm to 470nm)
%         idx = find(WL >= 420 & WL <= 470);
%         if ~isempty(idx)
%             sub_WL = WL(idx);
%             sub_Int = Intensity(idx);
%             center_wavelengths(k) = sum(sub_WL .* sub_Int) / sum(sub_Int);
%         else
%             center_wavelengths(k) = NaN;
%         end
% 
%         fprintf('%-40s | %-6.2f nm\n', filenames{k}, center_wavelengths(k));
% 
%         % Rescaling
%         max_val = max(Intensity);
%         if max_val == 0; max_val = 1; end 
%         all_normalized(:, k) = Intensity / max_val;
%         all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
%     end
% 
%     % --- Plotting ---
%     plot_spectra(all_WL, all_normalized, 'Normalized (Max=1)', labels);
%     plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);
% 
%     % 3. Optional CSV Outputs
%     answer = questdlg('Export processed data to two CSV files?', 'Export Data', 'Yes', 'No', 'No');
%     
%     if strcmp(answer, 'Yes')
%         [save_name, save_path] = uiputfile('Spectra_MaxUnity.csv', 'Save Normalized Data As');
%         if ~isequal(save_name, 0)
%             % File 1: Rescaled to Unity (Max=1)
%             export_csv_with_headers(fullfile(save_path, save_name), all_WL, all_normalized, filenames);
%             
%             % File 2: Rescaled to Common Integration Time
%             % Generate a default name for the second file
%             time_scaled_name = strrep(save_name, 'MaxUnity', 'TimeCorrected');
%             if strcmp(time_scaled_name, save_name) % if name didn't change, force a new one
%                 time_scaled_name = ['TimeCorrected_', save_name];
%             end
%             
%             [save_name2, save_path2] = uiputfile(time_scaled_name, 'Save Time-Corrected Data As');
%             if ~isequal(save_name2, 0)
%                 export_csv_with_headers(fullfile(save_path2, save_name2), all_WL, all_time_scaled, filenames);
%             end
%             
%             fprintf('\nBoth CSV files exported successfully.\n');
%         end
%     end
% end
% 
% % --- Helper Function: Export with Headers ---
% function export_csv_with_headers(full_file_path, wl, data_matrix, headers)
%     fid = fopen(full_file_path, 'w');
%     % Header line
%     fprintf(fid, 'Wavelength_nm');
%     for i = 1:length(headers)
%         fprintf(fid, ',%s', headers{i});
%     end
%     fprintf(fid, '\n');
%     fclose(fid);
%     
%     % Data block
%     export_matrix = [wl, data_matrix];
%     dlmwrite(full_file_path, export_matrix, '-append', 'delimiter', ',');
% end
% 
% % --- Helper Function: Plotting ---
% function plot_spectra(x, y, title_str, legs)
%     figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
%     plot(x, y, 'LineWidth', 2); 
%     grid on; xlim([440,460]);
%     xlabel('Wavelength (nm)', 'FontSize', 12);
%     ylabel('Intensity', 'FontSize', 12);
%     title(title_str, 'FontSize', 14);
%     legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
% end

% function batch_analyze_spectra_final()
%     % --- Setup Paths ---
%     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_16_LabData';
%     if ~exist(default_dir, 'dir'); default_dir = pwd; end
% 
%     % 1. Multi-select files
%     [filenames, start_path] = uigetfile({'*.txt', 'Avantes Spectra (*.txt)'}, ...
%         'Select Spectra Files', default_dir, 'MultiSelect', 'on');
%     
%     if isequal(filenames, 0); return; end
%     if ischar(filenames); filenames = {filenames}; end 
%     
%     num_files = length(filenames);
%     target_int_time = 100; 
%     
%     % Initialize containers
%     all_WL = [];
%     all_normalized = [];
%     all_time_scaled = [];
%     labels = cell(1, num_files);
%     center_wavelengths = zeros(1, num_files);
% 
%     fprintf('\n--- Processing Results ---\n');
%     fprintf('%-40s | %-20s\n', 'Filename', 'Center WL (420-470nm)');
%     fprintf('%s\n', repmat('-', 1, 65));
% 
%     for k = 1:num_files
%         full_path = fullfile(start_path, filenames{k});
%         labels{k} = strrep(filenames{k}, '_', ' ');
% 
%         % --- Extract Integration Time ---
%         int_time = 1; 
%         fid = fopen(full_path, 'r');
%         for i = 1:9
%             line_str = fgetl(fid);
%             if ~isempty(strfind(lower(line_str), 'integration time'))
%                 nums = sscanf(line_str, '%*[^0-9.]%f'); 
%                 if ~isempty(nums); int_time = nums(1); end
%             end
%         end
%         fclose(fid);
% 
%         % --- Read Numeric Data ---
%         data = dlmread(full_path, ';', 9, 0);
%         WL = data(:,1);
%         Intensity = data(:,2);
%         if k == 1; all_WL = WL; end 
% 
%         % 2. Calculate Center of Mass (420nm to 470nm)
%         idx = find(WL >= 420 & WL <= 470);
%         if ~isempty(idx)
%             sub_WL = WL(idx);
%             sub_Int = Intensity(idx);
%             % CoM Formula: sum(WL * Intensity) / sum(Intensity)
%             center_wavelengths(k) = sum(sub_WL .* sub_Int) / sum(sub_Int);
%         else
%             center_wavelengths(k) = NaN;
%         end
% 
%         % Print to terminal
%         fprintf('%-40s | %-6.2f nm\n', filenames{k}, center_wavelengths(k));
% 
%         % Rescaling
%         max_val = max(Intensity);
%         if max_val == 0; max_val = 1; end 
%         all_normalized(:, k) = Intensity / max_val;
%         all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
%     end
% 
%     % --- Plotting ---
%     plot_spectra(all_WL, all_normalized, 'Normalized (Max=1)', labels);
%     plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);
% 
%     % 3. Optional CSV Output
%     answer = questdlg('Would you like to generate a CSV output of the processed data?', ...
%         'Export Data', 'Yes', 'No', 'No');
%     
%     if strcmp(answer, 'Yes')
%         [save_name, save_path] = uiputfile('Processed_Spectra.csv', 'Save CSV As');
%         if ~isequal(save_name, 0)
%             % Combine WL, Normalized Data, and Time-Scaled Data
%             export_data = [all_WL, all_normalized, all_time_scaled];
%             csvwrite(fullfile(save_path, save_name), export_data);
%             fprintf('\nData exported to: %s\n', save_name);
%         end
%     end
% end
% 
% function plot_spectra(x, y, title_str, legs)
%     figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
%     % Increase LineWidth for visibility
%     p = plot(x, y, 'LineWidth', 2); 
%     grid on;
%     xlim([440,460]);
%     xlabel('Wavelength (nm)', 'FontSize', 12);
%     ylabel('Intensity', 'FontSize', 12);
%     title(title_str, 'FontSize', 14);
%     
%     % Move legend outside to the right
%     legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
% end


% function batch_analyze_Avantes_spectra()
%     % --- Setup Paths ---
%     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_16_LabData';
%     if ~exist(default_dir, 'dir'); default_dir = pwd; end
% 
%     % 1. Multi-select files
%     [filenames, start_path] = uigetfile({'*.txt', 'Avantes Spectra (*.txt)'}, ...
%         'Select Spectra Files', default_dir, 'MultiSelect', 'on');
%     
%     if isequal(filenames, 0); return; end
%     if ischar(filenames); filenames = {filenames}; end 
%     
%     num_files = length(filenames);
%     target_int_time = 100; 
%     
%     % Initialize containers
%     all_WL = [];
%     all_normalized = [];
%     all_time_scaled = [];
%     labels = cell(1, num_files);
% 
%     for k = 1:num_files
%         full_path = fullfile(start_path, filenames{k});
%         labels{k} = strrep(filenames{k}, '_', ' ');
% 
%         % --- Extract Integration Time from Header (Legacy Style) ---
%         int_time = 1; % Default
%         fid = fopen(full_path, 'r');
%         for i = 1:9
%             line_str = fgetl(fid);
%             if ~isempty(strfind(lower(line_str), 'integration time'))
%                 % Use sscanf to find the first number in the string
%                 nums = sscanf(line_str, '%*[^0-9.]%f'); 
%                 if ~isempty(nums); int_time = nums(1); end
%             end
%         end
%         fclose(fid);
% 
%         % --- Read Numeric Data (Using dlmread for older versions) ---
%         % dlmread(filename, delimiter, row_offset, column_offset)
%         data = dlmread(full_path, ';', 9, 0);
%         
%         WL = data(:,1);
%         Intensity = data(:,2);
% 
%         if k == 1; all_WL = WL; end 
% 
%         % 2. Rescale so max = 1
%         max_val = max(Intensity);
%         if max_val == 0; max_val = 1; end % Prevent divide by zero
%         all_normalized(:, k) = Intensity / max_val;
% 
%         % 4. Rescale to 100ms integration time
%         all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
%     end
% 
%     % --- Plotting ---
%     
%     % 3. Overlay Plot: Normalized
%     figure('Name', 'Normalized Spectra (Max=1)', 'Color', 'w');
%     plot(all_WL, all_normalized);
%     xlim([440,460]);
%     grid on; xlabel('Wavelength (nm)'); ylabel('Normalized Intensity');
%     title('Normalized Spectra: Peak Scaled to 1');
%     legend(labels, 'Interpreter', 'none');
% 
%     % 5. Overlay Plot: Rescaled to 100ms
%     figure('Name', 'Time-Corrected Spectra (100ms)', 'Color', 'w');
%     plot(all_WL, all_time_scaled);  
%     xlim([440,460]);
%     grid on; xlabel('Wavelength (nm)'); ylabel('Intensity (scaled to 100ms)');
%     title(sprintf('Spectra Scaled to %0.1f ms Integration Time', target_int_time));
%     legend(labels, 'Interpreter', 'none');
% 
%     fprintf('Done. Processed %d files using legacy readers.\n', num_files);
% end
%% Zemax POP Data Visualizer (Irradiance & Phase)
% Legacy & Modern MATLAB Compatible (R2018b and Later)

clear; clc; close all;

default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL\data_from_NECSEL';
if exist(default_path, 'dir')
    cd(default_path);
end

%% 1. Select the Zemax POP Text File
[filename, pathname] = uigetfile('*.txt', 'Select the Zemax POP Export File');
if isequal(filename, 0) || isequal(pathname, 0)
    disp('User canceled the operation.');
    return;
end
fullpath = fullfile(pathname, filename);

%% 2. Read File as Raw Binary and Convert to String
fid = fopen(fullpath, 'rb'); 
if fid == -1
    error('Could not open the selected file.');
end
rawBytes = fread(fid, [1, inf], 'uint8=>uint8');
fclose(fid);

% Convert bytes to character string. Handles Zemax Unicode (UTF-16LE) or ASCII/UTF-8
fileContent = native2unicode(rawBytes, 'UTF-16LE');
if isempty(fileContent) || ~contains(fileContent, 'Grid size')
    fileContent = native2unicode(rawBytes, 'UTF-8');
end

%% 3. Split into Lines and Parse Header
fileLines = textscan(fileContent, '%s', 'Delimiter', '\n', 'Whitespace', '');
fileLines = fileLines{1};

% Initialize variables to store header info with sensible defaults
titleStr = 'Zemax POP Data';
gridX = 512; gridY = 512; 
spacingX = 1; spacingY = 1;
wavelength = '';
peakValStr = 'N/A';
centerPhaseStr = 'N/A';
couplingEff = 'N/A';
isPhase = false; % Default assumption

% Parse through the first 20 rows of text using robust legacy scanning
numRowsToSearch = min(20, length(fileLines));
for row = 1:numRowsToSearch
    line = strtrim(fileLines{row});
    
    % 1. Explicitly Detect Data Type
    if contains(line, 'POP Irradiance Data') || contains(line, 'Total Irradiance surface')
        isPhase = false;
    elseif contains(line, 'POP Phase Data') || contains(line, 'Phase surface')
        isPhase = true;
    
    % 2. Extract Title
    elseif startsWith(line, 'Title:')
        titleStr = strtrim(strrep(line, 'Title:', ''));
    
    % 3. Extract Grid Size
    elseif contains(line, 'Grid size')
        % Look for the numbers directly after the text
        idx = strfind(line, ':');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%d by %d');
            if length(nums) == 2
                gridX = nums(1);
                gridY = nums(2);
            end
        end
        
    % 4. Extract Point Spacing
    elseif contains(line, 'Point spacing')
        idx = strfind(line, ':');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%f by %f');
            if length(nums) == 2
                spacingX = nums(1);
                spacingY = nums(2);
            end
        end
        
    % 5. Extract Wavelength
    elseif startsWith(line, 'Wavelength')
        wavelength = line;
        
    % 6. Extract Peak Irradiance (Only exists in Irradiance file)
    elseif contains(line, 'Peak Irradiance =')
        idx = strfind(line, '=');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%f');
            if ~isempty(nums)
                peakValStr = sprintf('%E', nums(1));
            end
        end
        
    % 7. Extract Center Phase (Only exists in Phase file)
    elseif contains(line, 'Center Phase =')
        idx = strfind(line, '=');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%f');
            if ~isempty(nums)
                centerPhaseStr = sprintf('%.4f', nums(1));
            end
        end
        
    % 8. Extract Fiber/Mode Coupling Efficiency Metrics
    elseif contains(line, 'Fiber Efficiency:')
        idx = strfind(line, 'Coupling');
        if ~isempty(idx)
            nums = sscanf(line(idx+8:end), '%f'); % Offset by length of 'Coupling'
            if ~isempty(nums)
                couplingEff = sprintf('%.2f%%', nums(1) * 100);
            end
        end
    end
end

% Set display strings based on detected file type
if isPhase
    typeStr = 'Phase Map';
    unitStr = 'Phase (radians)';
    subTitleStr = sprintf('Grid: %dx%d | Center Phase: %s rad | Zemax Coupling: %s', ...
                          gridX, gridY, centerPhaseStr, couplingEff);
else
    typeStr = 'Irradiance Map';
    unitStr = 'Irradiance (W/mm^2)';
    subTitleStr = sprintf('Grid: %dx%d | Peak: %s W/mm^2 | Zemax Coupling: %s', ...
                          gridX, gridY, peakValStr, couplingEff);
end

%% 4. Parse the Numeric Data Matrix
dataTextBlock = strjoin(fileLines(18:end), ' ');
allNumbers = sscanf(dataTextBlock, '%f');
expectedElements = gridX * gridY;

if isempty(allNumbers) || length(allNumbers) < expectedElements
    error('Could not parse data matrix. Expected %d elements, found %d.', ...
          expectedElements, length(allNumbers));
end

if length(allNumbers) > expectedElements
    allNumbers = allNumbers(1:expectedElements);
end

% Shape the flat array into the 2D grid (Zemax writes row-by-row)
popData = reshape(allNumbers, [gridX, gridY])';

%% 5. Generate Spatial Axes (Centered at 0,0)
x_axis = (((1:gridX) - (gridX + 1) / 2) * spacingX);
y_axis = (((1:gridY) - (gridY + 1) / 2) * spacingY);

%% 6. Plot the Heatmap with Context-Aware Styling
figure('Color', 'w', 'Position', [100, 100, 750, 600]);

imagesc(x_axis, y_axis, popData);
set(gca, 'YDir', 'normal'); 
axis image;                 
grid on;

if isPhase
    colormap(hsv); 
    cb = colorbar;
    ylabel(cb, unitStr, 'FontSize', 10);
    set(gca, 'GridColor', 'k', 'GridAlpha', 0.2); 
else
    colormap(jet); 
    cb = colorbar;
    ylabel(cb, unitStr, 'FontSize', 10);
    set(gca, 'GridColor', 'w', 'GridAlpha', 0.3); 
end

xlabel('X Position (mm)', 'FontSize', 11);
ylabel('Y Position (mm)', 'FontSize', 11);

titleLines = { ...
    sprintf('POP %s: %s', typeStr, titleStr), ...
    subTitleStr, ...
    wavelength ...
};
title(titleLines, 'Interpreter', 'none', 'FontSize', 12);

% 
% 
% %% Zemax POP Irradiance Map Visualizer
% % Compatible with MATLAB R2018b and later
% 
% clear; clc; close all;
% 
% default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL\';
% cd(default_path);
% 
% %% 1. Select the Zemax POP Text File
% [filename, pathname] = uigetfile('*.txt', 'Select the Zemax POP Irradiance File');
% if isequal(filename, 0) || isequal(pathname, 0)
%     disp('User canceled the operation.');
%     return;
% end
% fullpath = fullfile(pathname, filename);
% 
% %% 2. Read File as Raw Binary and Convert to String
% fid = fopen(fullpath, 'rb'); 
% if fid == -1
%     error('Could not open the selected file.');
% end
% rawBytes = fread(fid, [1, inf], 'uint8=>uint8');
% fclose(fid);
% 
% % Convert bytes to character string. 
% % If it's a Zemax Unicode export, 'UTF-16LE' handles it perfectly here.
% fileContent = native2unicode(rawBytes, 'UTF-16LE');
% 
% % Fallback: If it was already a normal ASCII/UTF-8 file, the line above might 
% % produce garbage text. Check if it looks right, otherwise interpret as UTF-8.
% if isempty(fileContent) || ~contains(fileContent, 'Grid size')
%     fileContent = native2unicode(rawBytes, 'UTF-8');
% end
% 
% %% 3. Split into Lines and Parse Header
% % Split the massive text block into individual cell rows by line breaks
% fileLines = textscan(fileContent, '%s', 'Delimiter', '\n', 'Whitespace', '');
% fileLines = fileLines{1};
% 
% % Initialize variables to store header info with sensible defaults
% titleStr = 'Zemax POP Data';
% gridX = 512; gridY = 512; 
% spacingX = 1; spacingY = 1;
% wavelength = '';
% peakIrrad = '';
% 
% % Parse through the first 17 rows of text
% numRowsToSearch = min(17, length(fileLines));
% for row = 1:numRowsToSearch
%     line = strtrim(fileLines{row});
%     
%     % Extract Title
%     if startsWith(line, 'Title:')
%         titleStr = strtrim(strrep(line, 'Title:', ''));
%     
%     % Extract Grid Size
%     elseif contains(line, 'Grid size (X by Y):')
%         tokens = regexp(line, 'Grid size \(X by Y\):\s*(\d+)\s*by\s*(\d+)', 'tokens');
%         if ~isempty(tokens)
%             gridX = str2double(tokens{1}{1});
%             gridY = str2double(tokens{1}{2});
%         end
%         
%     % Extract Point Spacing
%     elseif contains(line, 'Point spacing (X by Y):')
%         tokens = regexp(line, 'Point spacing \(X by Y\):\s*([\d\.E\+\-]+)\s*by\s*([\d\.E\+\-]+)', 'tokens');
%         if ~isempty(tokens)
%             spacingX = str2double(tokens{1}{1});
%             spacingY = str2double(tokens{1}{2});
%         end
%         
%     % Extract Wavelength
%     elseif startsWith(line, 'Wavelength')
%         wavelength = line;
%         
%     % Extract Peak Irradiance
%     elseif contains(line, 'Peak Irradiance =')
%         tokens = regexp(line, 'Peak Irradiance\s*=\s*([\d\.E\+\-]+)', 'tokens');
%         if ~isempty(tokens)
%             peakIrrad = tokens{1}{1};
%         end
%     end
% end
% 
% %% 4. Parse the Numeric Irradiance Data
% % Re-join the lines starting from line 18 back into a text block for fast scanning
% dataTextBlock = strjoin(fileLines(18:end), ' ');
% 
% % Scan all numeric values out of the block
% allNumbers = sscanf(dataTextBlock, '%f');
% 
% % Check if we got the expected amount of data
% expectedElements = gridX * gridY;
% 
% if isempty(allNumbers) || length(allNumbers) < expectedElements
%     error('Could not parse data matrix. Expected %d elements, found %d.', ...
%           expectedElements, length(allNumbers));
% end
% 
% % Truncate any trailing metadata if it exists
% if length(allNumbers) > expectedElements
%     allNumbers = allNumbers(1:expectedElements);
% end
% 
% % Shape the flat array into the 2D grid (Zemax writes row-by-row)
% irradianceData = reshape(allNumbers, [gridX, gridY])';
% 
% %% 5. Generate Spatial Axes (Centered at 0,0)
% % Zemax POP grids are centered on the optical axis
% x_axis = (((1:gridX) - (gridX + 1) / 2) * spacingX);
% y_axis = (((1:gridY) - (gridY + 1) / 2) * spacingY);
% 
% %% 6. Plot the Heatmap
% figure('Color', 'w', 'Position', [100, 100, 750, 600]);
% 
% % Heatmap presentation
% imagesc(x_axis, y_axis, irradianceData);
% set(gca, 'YDir', 'normal'); % Flips Y-axis so standard Cartesian math applies
% 
% % Visual Styling
% colormap(jet); 
% cb = colorbar;
% ylabel(cb, 'Irradiance (W/mm^2)', 'FontSize', 10);
% axis image; % Fixes the 1:1 aspect ratio aspect of the beam shape
% 
% % Labels and Dynamic Title
% xlabel('X Position (mm)', 'FontSize', 11);
% ylabel('Y Position (mm)', 'FontSize', 11);
% 
% titleLines = { ...
%     sprintf('POP Irradiance Map: %s', titleStr), ...
%     sprintf('Grid: %dx%d | Peak: %s W/mm^2', gridX, gridY, peakIrrad), ...
%     wavelength ...
% };
% title(titleLines, 'Interpreter', 'none', 'FontSize', 12);
% 
% grid on;
% set(gca, 'GridColor', 'w', 'GridAlpha', 0.3);

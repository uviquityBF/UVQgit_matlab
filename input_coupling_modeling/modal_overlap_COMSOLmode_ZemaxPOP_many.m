%% Batch Waveguide Coupling & Mode Overlap Calculator
% Compatible with MATLAB R2018b and later

clear; clc; close all;

% Setup default directory if it exists
default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL';
if exist(default_path, 'dir'), cd(default_path); end

%% ==========================================
%% 1. LOAD AND PARSE STATIC ZEMAX POP DATA
%% ==========================================
fprintf('--- Step 1: Loading Static Zemax POP Data ---\n');

% 1a. Load Zemax Irradiance File
[fileI, pathI] = uigetfile('*.txt', 'Select the Zemax POP IRRADIANCE File');
if isequal(fileI,0), disp('Operation canceled.'); return; end
fid = fopen(fullfile(pathI, fileI), 'rb');
rawBytes = fread(fid, [1, inf], 'uint8=>uint8'); fclose(fid);
txtI = native2unicode(rawBytes, 'UTF-16LE');
if isempty(txtI) || ~contains(txtI, 'Grid size'), txtI = native2unicode(rawBytes, 'UTF-8'); end
linesI = textscan(txtI, '%s', 'Delimiter', '\n', 'Whitespace', ''); linesI = linesI{1};

% 1b. Load Zemax Phase File
[fileP, pathP] = uigetfile('*.txt', 'Select the matching Zemax POP PHASE File');
if isequal(fileP,0), disp('Operation canceled.'); return; end
fid = fopen(fullfile(pathP, fileP), 'rb');
rawBytes = fread(fid, [1, inf], 'uint8=>uint8'); fclose(fid);
txtP = native2unicode(rawBytes, 'UTF-16LE');
if isempty(txtP) || ~contains(txtP, 'Grid size'), txtP = native2unicode(rawBytes, 'UTF-8'); end
linesP = textscan(txtP, '%s', 'Delimiter', '\n', 'Whitespace', ''); linesP = linesP{1};

% Parse Grid Size & Spacing from Irradiance Header
zGridX = 512; zGridY = 512; zSpacingX = 1; zSpacingY = 1;
for r = 1:min(20, length(linesI))
    line = strtrim(linesI{r});
    if contains(line, 'Grid size')
        idx = strfind(line, ':');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%d by %d');
            zGridX = nums(1); zGridY = nums(2);
        end
    elseif contains(line, 'Point spacing')
        idx = strfind(line, ':');
        if ~isempty(idx)
            nums = sscanf(line(idx+1:end), '%f by %f');
            zSpacingX = nums(1); zSpacingY = nums(2);
        end
    end
end

% Parse Numeric Matrices
numI = sscanf(strjoin(linesI(18:end), ' '), '%f');
numP = sscanf(strjoin(linesP(18:end), ' '), '%f');
I_zemax_raw = reshape(numI(1:(zGridX*zGridY)), [zGridX, zGridY])';
Phi_zemax_raw = reshape(numP(1:(zGridX*zGridY)), [zGridX, zGridY])';

% Generate Zemax Axes (Centered at 0,0 in mm)
zX = (((1:zGridX) - (zGridX + 1) / 2) * zSpacingX);
zY = (((1:zGridY) - (zGridY + 1) / 2) * zSpacingY);
[Z_Xmesh, Z_Ymesh] = meshgrid(zX, zY);

%% ==========================================
%% 2. SELECT MULTIPLE COMSOL WAVEGUIDE FILES
%% ==========================================
fprintf('\n--- Step 2: Selecting Waveguide Mode Files ---\n');

[comsolFiles, pathC] = uigetfile({'*.txt;*.dat', 'COMSOL Files (*.txt, *.dat)'}, ...
                                  'Select ONE or MULTIPLE COMSOL Mode Files', ...
                                  'MultiSelect', 'on');

if isequal(comsolFiles, 0)
    disp('Operation canceled.');
    return;
end

% Standardize input format into a cell array even if only one file is selected
if ~iscell(comsolFiles)
    comsolFiles = {comsolFiles};
end

numFiles = length(comsolFiles);
fprintf('Found %d files queued for batch evaluation.\n', numFiles);

% Initialize storage structures for our final reporting table
results_FileName = cell(numFiles, 1);
results_Coupling = zeros(numFiles, 1);

% Configure text import properties once to speed up processing
opts = delimitedTextImportOptions('NumVariables', 4, 'Delimiter', ',', 'CommentStyle', '%');
opts.VariableTypes = {'double', 'double', 'double', 'double'};

%% ==========================================
%% 3. BATCH PROCESSING LOOP
%% ==========================================
for fIdx = 1:numFiles
    currentFile = comsolFiles{fIdx};
    fprintf('\nProcessing File %d/%d: %s\n', fIdx, numFiles, currentFile);
    
    % 3a. Read and clean current COMSOL profile
    dataTable = readtable(fullfile(pathC, currentFile), opts);
    cX_mm = dataTable.Var1 / 1000;   % µm to mm
    cY_mm = dataTable.Var2 / 1000;   % µm to mm
    cNormE = dataTable.Var3;
    cPowerZ = dataTable.Var4;
    
    badRows = isnan(cX_mm) | isinf(cX_mm) | isnan(cY_mm) | isinf(cY_mm) | ...
               isnan(cNormE) | isinf(cNormE) | isnan(cPowerZ) | isinf(cPowerZ);
    cX_mm(badRows) = []; cY_mm(badRows) = []; cNormE(badRows) = []; cPowerZ(badRows) = [];
    
    % 3b. Determine 1% Intensity Bounding Box for this file combination
    maxI = max(I_zemax_raw(:));
    maskI = I_zemax_raw >= (0.01 * maxI);
    zX_active = Z_Xmesh(maskI); zY_active = Z_Ymesh(maskI);
    
    maxP = max(cPowerZ);
    maskC = cPowerZ >= (0.01 * maxP);
    cX_active = cX_mm(maskC); cY_active = cY_mm(maskC);
    
    minX = min([min(zX_active), min(cX_active)]);
    maxX = max([max(zX_active), max(cX_active)]);
    minY = min([min(zY_active), min(cY_active)]);
    maxY = max([max(zY_active), max(cY_active)]);
    
    % Generate local high-density uniform grid points
    gridPoints = 512;
    intX = linspace(minX, maxX, gridPoints);
    intY = linspace(minY, maxY, gridPoints);
    [Int_Xmesh, Int_Ymesh] = meshgrid(intX, intY);
    
    % 3c. Interpolation and Field Reconstruction
    I_int = griddata(zX, zY, I_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
    Phi_int = griddata(zX, zY, Phi_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
    I_int(isnan(I_int)) = 0; Phi_int(isnan(Phi_int)) = 0;
    
    E_laser = sqrt(I_int) .* exp(1i * Phi_int);
    
    E_chip = griddata(cX_mm, cY_mm, cNormE, Int_Xmesh, Int_Ymesh, 'linear');
    E_chip(isnan(E_chip)) = 0;
    
    % 3d. Overlap Calculation
    numerator = abs(trapz(intY, trapz(intX, E_laser .* conj(E_chip), 2)))^2;
    denom_laser = trapz(intY, trapz(intX, abs(E_laser).^2, 2));
    denom_chip  = trapz(intY, trapz(intX, abs(E_chip).^2, 2));
    
    coupling_efficiency = numerator / (denom_laser * denom_chip);
    
    % Save data point entries for the table output
    results_FileName{fIdx} = currentFile;
    results_Coupling(fIdx) = coupling_efficiency * 100; % Percentage
    
    % 3e. Generate and Save Visualizations
    % Generates figures sequentially. Figures are explicitly titled to avoid rewriting window overlaps.
    figHandle = figure('Color', 'w', 'Position', [50 + (fIdx*15), 50 + (fIdx*15), 1400, 450]);
    
    % Left Plot: Laser Profile
    subplot(1, 3, 1);
    imagesc(intX*1000, intY*1000, abs(E_laser));
    set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
    xlabel('X Position (µm)'); ylabel('Y Position (µm)');
    title('Input Laser Field |E|');
    
    % Middle Plot: Current Waveguide Target
    subplot(1, 3, 2);
    imagesc(intX*1000, intY*1000, abs(E_chip));
    set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
    xlabel('X Position (µm)'); ylabel('Y Position (µm)');
    title('Waveguide Target Field |E|');
    
    % Right Plot: Laser Wave Phase
    subplot(1, 3, 3);
    imagesc(intX*1000, intY*1000, angle(E_laser));
    set(gca, 'YDir', 'normal'); axis image; colormap(hsv); colorbar; grid on;
    xlabel('X Position (µm)'); ylabel('Y Position (µm)');
    title('Input Laser Phase (rad)');
    
    % Descriptive main visual text header
    annotation('textbox', [0, 0.91, 1, 0.09], ...
               'String', sprintf('File [%d/%d]: %s\\nCalculated Coupling: %.2f%%', fIdx, numFiles, currentFile, results_Coupling(fIdx)), ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
               'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
           
    drawnow; % Forces MATLAB to refresh UI and draw the figure immediately
end

%% ==========================================
%% 4. COMPILE, DISPLAY, AND EXPORT SUMMARY TABLE
%% ==========================================
fprintf('\n--- Step 4: Displaying Final Results Summary Table ---\n');

% Construct the final MATLAB Dataset Table
SummaryTable = table(results_FileName, results_Coupling, ...
                     'VariableNames', {'COMSOL_File_Name', 'Overlap_Efficiency_Percent'});

% Sort results automatically by performance (highest coupling first)
SummaryTable = sortrows(SummaryTable, 'Overlap_Efficiency_Percent', 'descend');

% Print results table explicitly to the command window layout
disp(' ');
disp('================================================================');
disp('                FINAL COUPLING EFFICIENCY SUMMARY TABLE         ');
disp('================================================================');
disp(SummaryTable);
disp('================================================================');

% Optional File Export: Uncomment line below to write directly to a local CSV file sheet
% writetable(SummaryTable, 'Batch_Coupling_Analysis_Output.csv');
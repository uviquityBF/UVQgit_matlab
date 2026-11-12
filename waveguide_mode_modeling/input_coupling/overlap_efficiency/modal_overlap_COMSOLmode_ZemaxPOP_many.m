%% Batch Waveguide Coupling & Mode Overlap Calculator
% Compatible with MATLAB R2018b and later

addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\waveguide');
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
POP_I = load_zemax_pop_txt(fullfile(pathI, fileI));

% 1b. Load Zemax Phase File
[fileP, pathP] = uigetfile('*.txt', 'Select the matching Zemax POP PHASE File');
if isequal(fileP,0), disp('Operation canceled.'); return; end
POP_P = load_zemax_pop_txt(fullfile(pathP, fileP));

% Grid Size & Spacing taken from the Irradiance file header
zGridX = POP_I.gridX; zGridY = POP_I.gridY;
I_zemax_raw = POP_I.data;
Phi_zemax_raw = POP_P.data;

% Zemax Axes (Centered at 0,0 in mm)
zX = POP_I.x_axis;
zY = POP_I.y_axis;
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

%% ==========================================
%% 3. BATCH PROCESSING LOOP
%% ==========================================
for fIdx = 1:numFiles
    currentFile = comsolFiles{fIdx};
    fprintf('\nProcessing File %d/%d: %s\n', fIdx, numFiles, currentFile);

    % 3a. Read and clean current COMSOL profile
    [cX_um, cY_um, cNormE, cPowerZ] = load_comsol_mode_profile(fullfile(pathC, currentFile));
    cX_mm = cX_um / 1000;   % um to mm
    cY_mm = cY_um / 1000;   % um to mm

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
    xlabel('X Position (�m)'); ylabel('Y Position (�m)');
    title('Input Laser Field |E|');
    
    % Middle Plot: Current Waveguide Target
    subplot(1, 3, 2);
    imagesc(intX*1000, intY*1000, abs(E_chip));
    set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
    xlabel('X Position (�m)'); ylabel('Y Position (�m)');
    title('Waveguide Target Field |E|');
    
    % Right Plot: Laser Wave Phase
    subplot(1, 3, 3);
    imagesc(intX*1000, intY*1000, angle(E_laser));
    set(gca, 'YDir', 'normal'); axis image; colormap(hsv); colorbar; grid on;
    xlabel('X Position (�m)'); ylabel('Y Position (�m)');
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
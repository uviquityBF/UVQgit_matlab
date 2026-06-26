%% Master Waveguide Coupling & Mode Overlap Calculator
% Compatible with MATLAB R2018b and later

clear; clc; close all;

% Setup default directory if it exists
default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL';
if exist(default_path, 'dir'), cd(default_path); end

%% ==========================================
%% 1. LOAD AND PARSE ZEMAX POP DATA (LASER INPUT)
%% ==========================================
fprintf('--- Step 1: Loading Zemax POP Data ---\n');

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

% Parse Grid Size & Spacing from Irradiance Header (Assuming Phase matches)
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
%% 2. LOAD AND CLEAN COMSOL DATA (CHIP MODE)
%% ==========================================
fprintf('\n--- Step 2: Loading COMSOL Waveguide Data ---\n');

[fileC, pathC] = uigetfile({'*.txt;*.dat'}, 'Select the COMSOL Waveguide Export File');
if isequal(fileC,0), disp('Operation canceled.'); return; end

opts = delimitedTextImportOptions('NumVariables', 4, 'Delimiter', ',', 'CommentStyle', '%');
opts.VariableTypes = {'double', 'double', 'double', 'double'};
dataTable = readtable(fullfile(pathC, fileC), opts);

cX_mm = dataTable.Var1 / 1000;   % Convert µm to mm to align with Zemax
cY_mm = dataTable.Var2 / 1000;   % Convert µm to mm to align with Zemax
cNormE = dataTable.Var3;
cPowerZ = dataTable.Var4;

% Data Cleaning (Purge NaNs/Infs)
badRows = isnan(cX_mm) | isinf(cX_mm) | isnan(cY_mm) | isinf(cY_mm) | ...
           isnan(cNormE) | isinf(cNormE) | isnan(cPowerZ) | isinf(cPowerZ);
cX_mm(badRows) = []; cY_mm(badRows) = []; cNormE(badRows) = []; cPowerZ(badRows) = [];

%% ==========================================
%% 3. DEFINE DYNAMIC DENSE INTEGRATION GRID (1% Threshold)
%% ==========================================
fprintf('\n--- Step 3: Determining 1%% Intensity Bounding Box ---\n');

% Find spatial limits where Zemax Intensity is > 1% of peak
maxI = max(I_zemax_raw(:));
maskI = I_zemax_raw >= (0.01 * maxI);
zX_active = Z_Xmesh(maskI);
zY_active = Z_Ymesh(maskI);

% Find spatial limits where COMSOL Power Intensity is > 1% of peak
maxP = max(cPowerZ);
maskC = cPowerZ >= (0.01 * maxP);
cX_active = cX_mm(maskC);
cY_active = cY_mm(maskC);

% Define the absolute outer envelope boundary encompassing both active domains
minX = min([min(zX_active), min(cX_active)]);
maxX = max([max(zX_active), max(cX_active)]);
minY = min([min(zY_active), min(cY_active)]);
maxY = max([max(zY_active), max(cY_active)]);

fprintf('Integration Domain Restricted To:\n');
fprintf('  X: [%.4f to %.4f] mm\n', minX, maxX);
fprintf('  Y: [%.4f to %.4f] mm\n', minY, maxY);

% Generate the pristine high-density uniform Integration Grid over this limited window
gridPoints = 512;
intX = linspace(minX, maxX, gridPoints);
intY = linspace(minY, maxY, gridPoints);
[Int_Xmesh, Int_Ymesh] = meshgrid(intX, intY);

%% ==========================================
%% 4. INTERPOLATION AND COMPLEX FIELD RECONSTRUCTION
%% ==========================================
fprintf('\n--- Step 4: Reconstructing Fields on Dense Integration Grid ---\n');

% Interpolate Zemax Raw Intensity and Phase onto the localized grid
I_int = griddata(zX, zY, I_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
Phi_int = griddata(zX, zY, Phi_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
I_int(isnan(I_int)) = 0; Phi_int(isnan(Phi_int)) = 0;

% Construct the Complex Electric Field Vector from the Laser Input
E_laser = sqrt(I_int) .* exp(1i * Phi_int);

% Interpolate COMSOL Scattered Mode Profile onto the exact same grid
E_chip = griddata(cX_mm, cY_mm, cNormE, Int_Xmesh, Int_Ymesh, 'linear');
E_chip(isnan(E_chip)) = 0;

%% ==========================================
%% 5. CALCULATE MODE COUPLING EFFICIENCY (OVERLAP INTEGRAL)
%% ==========================================
fprintf('\n--- Step 5: Calculating Overlap Integral ---\n');

% Numerator: |Double Integral( E_laser * conj(E_chip) dx dy) |^2
% Since grid spacing is uniform, dx*dy constants cancel out in the division.
numerator = abs(trapz(intY, trapz(intX, E_laser .* conj(E_chip), 2)))^2;

% Denominator: Integral(|E_laser|^2) * Integral(|E_chip|^2)
denom_laser = trapz(intY, trapz(intX, abs(E_laser).^2, 2));
denom_chip  = trapz(intY, trapz(intX, abs(E_chip).^2, 2));

coupling_efficiency = numerator / (denom_laser * denom_chip);

fprintf('\n======================================\n');
fprintf('  SUCCESS: OVERLAP COUPLING EFFICIENCY\n');
fprintf('  Calculated Coupling: %.2f%%\n', coupling_efficiency * 100);
fprintf('======================================\n');

%% ==========================================
%% 6. PLOT COMPARATIVE RESULTS
%% ==========================================
figure('Color', 'w', 'Position', [50, 50, 1400, 450]);

% Left Plot: Laser Field Magnitude on Integration Grid
subplot(1, 3, 1);
imagesc(intX*1000, intY*1000, abs(E_laser));
set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
xlabel('X Position (µm)'); ylabel('Y Position (µm)');
title('Input Laser Field |E| (Interpolated)');

% Middle Plot: Waveguide Field Magnitude on Integration Grid
subplot(1, 3, 2);
imagesc(intX*1000, intY*1000, abs(E_chip));
set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
xlabel('X Position (µm)'); ylabel('Y Position (µm)');
title('Waveguide Target Field |E| (Interpolated)');

% Right Plot: Laser Field Phase on Integration Grid
subplot(1, 3, 3);
imagesc(intX*1000, intY*1000, angle(E_laser));
set(gca, 'YDir', 'normal'); axis image; colormap(hsv); colorbar; grid on;
xlabel('X Position (µm)'); ylabel('Y Position (µm)');
title('Input Laser Phase Angle (rad)');

annotation('textbox', [0, 0.91, 1, 0.09], ...
           'String', sprintf('Master Overlap Efficiency Analysis\\nCalculated Coupling: %.2f%%', coupling_efficiency * 100), ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
           'FontSize', 14, 'FontWeight', 'bold');
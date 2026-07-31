function POP = load_zemax_pop_txt(filepath)
% Loads a Zemax POP text export (irradiance or phase). Handles Zemax's
% UTF-16LE encoding with a UTF-8 fallback, scans the header for Grid
% size / Point spacing / Wavelength / Peak Irradiance / Center Phase /
% Fiber Efficiency, auto-detects whether the file holds irradiance or
% phase data, and reshapes the numeric block into a 2D grid with
% centered spatial axes.
%
% USAGE:
%   POP = load_zemax_pop_txt(filepath);
%
% Returns a struct:
%   POP.data            [gridY x gridX] numeric grid (irradiance or phase)
%   POP.isPhase          true if this file was detected as phase data
%   POP.gridX, gridY     grid dimensions
%   POP.spacingX, spacingY
%   POP.x_axis, y_axis   centered spatial axes (row vectors)
%   POP.title            Zemax "Title:" string
%   POP.wavelength        raw "Wavelength..." header line
%   POP.peakValStr        "Peak Irradiance" string (Irradiance files only)
%   POP.centerPhaseStr    "Center Phase" string (Phase files only)
%   POP.couplingEff       "Fiber Efficiency" coupling string, if present

    fid = fopen(filepath, 'rb');
    if fid == -1
        error('Could not open the selected file: %s', filepath);
    end
    rawBytes = fread(fid, [1, inf], 'uint8=>uint8');
    fclose(fid);

    fileContent = native2unicode(rawBytes, 'UTF-16LE');
    if isempty(fileContent) || ~contains(fileContent, 'Grid size')
        fileContent = native2unicode(rawBytes, 'UTF-8');
    end

    fileLines = textscan(fileContent, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fileLines = fileLines{1};

    POP.title = 'Zemax POP Data';
    POP.gridX = 512; POP.gridY = 512;
    POP.spacingX = 1; POP.spacingY = 1;
    POP.wavelength = '';
    POP.peakValStr = 'N/A';
    POP.centerPhaseStr = 'N/A';
    POP.couplingEff = 'N/A';
    POP.isPhase = false;

    numRowsToSearch = min(20, length(fileLines));
    for row = 1:numRowsToSearch
        line = strtrim(fileLines{row});

        if contains(line, 'POP Irradiance Data') || contains(line, 'Total Irradiance surface')
            POP.isPhase = false;
        elseif contains(line, 'POP Phase Data') || contains(line, 'Phase surface')
            POP.isPhase = true;
        elseif startsWith(line, 'Title:')
            POP.title = strtrim(strrep(line, 'Title:', ''));
        elseif contains(line, 'Grid size')
            idx = strfind(line, ':');
            if ~isempty(idx)
                nums = sscanf(line(idx+1:end), '%d by %d');
                if length(nums) == 2
                    POP.gridX = nums(1);
                    POP.gridY = nums(2);
                end
            end
        elseif contains(line, 'Point spacing')
            idx = strfind(line, ':');
            if ~isempty(idx)
                nums = sscanf(line(idx+1:end), '%f by %f');
                if length(nums) == 2
                    POP.spacingX = nums(1);
                    POP.spacingY = nums(2);
                end
            end
        elseif startsWith(line, 'Wavelength')
            POP.wavelength = line;
        elseif contains(line, 'Peak Irradiance =')
            idx = strfind(line, '=');
            if ~isempty(idx)
                nums = sscanf(line(idx+1:end), '%f');
                if ~isempty(nums)
                    POP.peakValStr = sprintf('%E', nums(1));
                end
            end
        elseif contains(line, 'Center Phase =')
            idx = strfind(line, '=');
            if ~isempty(idx)
                nums = sscanf(line(idx+1:end), '%f');
                if ~isempty(nums)
                    POP.centerPhaseStr = sprintf('%.4f', nums(1));
                end
            end
        elseif contains(line, 'Fiber Efficiency:')
            idx = strfind(line, 'Coupling');
            if ~isempty(idx)
                nums = sscanf(line(idx+8:end), '%f'); % Offset by length of 'Coupling'
                if ~isempty(nums)
                    POP.couplingEff = sprintf('%.2f%%', nums(1) * 100);
                end
            end
        end
    end

    dataTextBlock = strjoin(fileLines(18:end), ' ');
    allNumbers = sscanf(dataTextBlock, '%f');
    expectedElements = POP.gridX * POP.gridY;

    if isempty(allNumbers) || length(allNumbers) < expectedElements
        error('Could not parse data matrix. Expected %d elements, found %d.', ...
              expectedElements, length(allNumbers));
    end
    if length(allNumbers) > expectedElements
        allNumbers = allNumbers(1:expectedElements);
    end

    % Zemax writes row-by-row.
    POP.data = reshape(allNumbers, [POP.gridX, POP.gridY])';

    POP.x_axis = (((1:POP.gridX) - (POP.gridX + 1) / 2) * POP.spacingX);
    POP.y_axis = (((1:POP.gridY) - (POP.gridY + 1) / 2) * POP.spacingY);
end

function DAT = read_spectrum_file(filepath, format)
% Reads a spectrometer text-file export into a common struct, covering the
% two spectrometer text formats used across this repo:
%
%   'oceanoptics' - Ocean Optics / QEPro / Maya (SpectraSuite/OceanView).
%                   Header block of metadata lines ending in a line
%                   containing "Begin Spectral Data", then two
%                   whitespace-delimited columns (wavelength, counts).
%                   Integration time is on a line containing "Integration
%                   Time" (value in seconds); "Scans to average" gives the
%                   scan count if present.
%
%   'avantes'     - Avantes ASCII export. Fixed 9-line header, then
%                   ';'-delimited columns (wavelength, intensity, dark,
%                   reference). Integration time is on a header line
%                   containing "integration time" (case-insensitive),
%                   pulled out via a trailing-number regex.
%
% USAGE:
%   DAT = read_spectrum_file(filepath, 'oceanoptics');
%   DAT = read_spectrum_file(filepath, 'avantes');
%   DAT = read_spectrum_file(filepath);            % auto-detect format
%
% DAT fields (unused ones for a given format are left as NaN/[]):
%   DAT.wl        - wavelength [nm]
%   DAT.counts    - counts (oceanoptics) or intensity (avantes)
%   DAT.dark      - dark column (avantes only)
%   DAT.reference - reference column (avantes only)
%   DAT.Tint      - integration time, native units from the file
%                   (seconds for oceanoptics, whatever the avantes header
%                   states -- callers have historically treated it as ms)
%   DAT.Nscans    - "Scans to average" (oceanoptics only, NaN if absent)
%   DAT.format    - the format actually used ('oceanoptics'/'avantes')
%
% This consolidates the file-reading logic previously duplicated (with
% minor variations/bugs) across QEPro_Spectra_Analysis.m,
% QEPro_Spectra_Animation.m, QEPro_OpenMany_Spectra.m,
% batch_analyze_spectra_multi_type.m, batch_open_avantes_spectra.m,
% spectrometer_CAL.m, APPLY_CAL.m, and APPLY_CAL_many.m.

    if nargin < 2 || isempty(format)
        format = detect_spectrum_format(filepath);
    end

    switch lower(format)
        case {'oceanoptics', 'qepro', 'maya'}
            DAT = read_oceanoptics(filepath);
            DAT.format = 'oceanoptics';
        case 'avantes'
            DAT = read_avantes(filepath);
            DAT.format = 'avantes';
        otherwise
            error('read_spectrum_file:unknownFormat', ...
                'Unknown format "%s" (expected ''oceanoptics'' or ''avantes'').', format);
    end
end


function format = detect_spectrum_format(filepath)
% Scans the first few lines for a reliable signal: Avantes headers are
% ';'-delimited, Ocean Optics headers say so in plain text.
    fid = fopen(filepath, 'r');
    if fid < 0
        error('read_spectrum_file:cannotOpen', 'Cannot open file: %s', filepath);
    end

    format = 'oceanoptics'; % default
    for k = 1:15
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if numel(strfind(line, ';')) >= 2
            format = 'avantes';
            break;
        end
        if contains(line, 'Begin Spectral Data') || contains(line, 'Integration Time')
            format = 'oceanoptics';
            break;
        end
    end
    fclose(fid);
end


function DAT = read_oceanoptics(filepath)
    fid = fopen(filepath, 'r');
    if fid < 0
        error('read_spectrum_file:cannotOpen', 'Cannot open file: %s', filepath);
    end

    DAT.Tint = NaN;
    DAT.Nscans = NaN;
    max_header_lines = 60; % real QEPro/Maya headers run ~13-14 lines
    found_marker = false;
    for k = 1:max_header_lines
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if contains(line, 'Integration Time') || contains(line, '(sec):')
            tok = regexp(line, '([\d.]+[Ee]?[+-]?\d*)\s*$', 'match', 'once');
            if ~isempty(tok), DAT.Tint = str2double(tok); end
        elseif contains(line, 'Scans to average')
            tok = regexp(line, '(\d+)\s*$', 'match', 'once');
            if ~isempty(tok), DAT.Nscans = str2double(tok); end
        elseif contains(line, 'Begin Spectral Data')
            found_marker = true;
            break;
        end
    end

    if ~found_marker
        % Fallback for headers without the literal marker: rewind and
        % scan for the first line that looks like numeric data.
        frewind(fid);
        for k = 1:max_header_lines
            pos = ftell(fid);
            line = fgetl(fid);
            if ~ischar(line)
                break;
            end
            if ~isempty(line) && (isstrprop(line(1), 'digit') || line(1) == '-')
                fseek(fid, pos, 'bof');
                break;
            end
        end
    end

    raw = textscan(fid, '%f %f');
    fclose(fid);

    DAT.wl = raw{1};
    DAT.counts = raw{2};
    DAT.dark = [];
    DAT.reference = [];
end


function DAT = read_avantes(filepath)
    n_header = 9;

    fid = fopen(filepath, 'r');
    if fid < 0
        error('read_spectrum_file:cannotOpen', 'Cannot open file: %s', filepath);
    end

    DAT.Tint = NaN;
    for k = 1:n_header
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if contains(lower(line), 'integration time')
            nums = sscanf(line, '%*[^0-9.]%f');
            if ~isempty(nums), DAT.Tint = nums(1); end
        end
    end
    fclose(fid);

    data = dlmread(filepath, ';', n_header, 0);
    DAT.wl = data(:, 1);
    DAT.counts = data(:, 2);
    if size(data, 2) >= 3
        DAT.dark = data(:, 3);
    else
        DAT.dark = [];
    end
    if size(data, 2) >= 4
        DAT.reference = data(:, 4);
    else
        DAT.reference = [];
    end
    DAT.Nscans = NaN;
end

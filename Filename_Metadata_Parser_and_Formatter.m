%% Filename Metadata Parser and Formatter
clear; clc;

% --- Configuration ---
inputDir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_23_LabData\SHG Search 50um\runG_am156_sd57_id98\';  % Folder containing source files
outputDir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_23_LabData\SHG Search 50um\copy_temp'; % Where the new files will go
inputExt  = '.txt';

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

files = dir(fullfile(inputDir, ['*', inputExt]));

for i = 1:length(files)
    oldName0 = files(i).name;
    oldName = oldName0;
    oldName = excise_Tag(oldName ,'RunG_am156_sd57');
    oldName = excise_Tag(oldName ,'slit=50um');
        
    % 1. Robust Date Extraction using Regex
    % Pattern looks for: 
    % ^\d{6} (start of string, 6 digits) OR 
    % ^\d{4}_\d{2}_\d{2} (start of string, YYYY_MM_DD)
    dateMatch = regexp(oldName, '^(\d{6}|\d{4}_\d{2}_\d{2})', 'match', 'once');
    
    if isempty(dateMatch)
        fprintf('Skipping %s: No valid date found at start.\n', oldName);
    end
    
    datePart = dateMatch;
    
    % 2. Get the "ID" part (the part immediately following the date)
    % We remove the date and the first underscore from the original name
    remainderAfterDate = extractAfter(oldName, [datePart, '_']);
  
    % Split the remainder to get the WG/ID part
    parts = split(remainderAfterDate, '_');
    wgPart = parts{1}; % This will be 'WG118'
    
    % --- Your existing ID conversion logic ---
    idPart = strrep(wgPart, 'WG', 'id');
    NumPerGroup = 3;
    tmp = split(idPart,'id');
    idnum = str2double(tmp{2}); % str2double is slightly safer than str2num
    idnum_new = (floor((idnum-1)/NumPerGroup)+1) + 0.1*(mod(idnum-1,NumPerGroup)+1);
    idPart_new = ['id', num2str(idnum_new)];
    idPart = idPart_new;
    
    % --- Reconstruct the start ---
    newNameStart = sprintf('%s_RunG_%s', datePart, idPart);
    
    % Calculate where the "remainingName" starts (everything after WGPart)
    % We search for the first underscore after the WG string in the remainder
    remainingName = extractAfter(remainderAfterDate, [wgPart, '_']);
    if isempty(remainingName), remainingName = ''; end % Handle short names
    % 3. Reconstruct the filename
    % We insert 'RunG' at index 2, then the modified ID part
    
    % Construct the beginning: Date_RunG_id118
%     newNameStart = sprintf('%s_RunG_%s', datePart, idPart);
    newNameStart = sprintf('%s_%s', datePart, idPart);
    
    % Get the rest of the original name after the "WG" part
    % We find where the second underscore was and take everything after it
    firstTwoPartsLength = length(datePart) + length(wgPart) + 2; % +2 for the underscores
    remainingName = oldName(firstTwoPartsLength:end);
    
    %modifications 
    remainingName = strrep(remainingName,'TiSWL=', 'pump');
    remainingName = strrep(remainingName,'1000ms', '1s');    
    newName = [newNameStart, remainingName];
    
    
    % 4. File Operations
    fullInputPath  = fullfile(inputDir, oldName0);
    fullOutputPath = fullfile(outputDir, newName);
    
    % Using 'copyfile' instead of read/write since you mentioned "copy to a format"
    % This preserves the internal data exactly as is.
    [status, msg] = copyfile(fullInputPath, fullOutputPath);
    
    if status
        fprintf('Success: %s \n      -> %s\n', oldName, newName);
    else
        fprintf('Error copying %s: %s\n', oldName, msg);
    end
end

% --- Metadata Tag Removal Section ---
function newstring = excise_Tag(oldstring,strToRemove)
tagToRemove = strToRemove;

    if contains(oldstring, tagToRemove)
        newstring = oldstring;
        % Option A: Remove the tag and the trailing underscore
        newstring = strrep(newstring, [tagToRemove, '_'], '');

        % Option B: If the tag is at the very end, the above might fail, 
        % so we do a secondary cleanup just in case:
        newstring = strrep(newstring, tagToRemove, '');

        % Clean up any accidental triple underscores created during the swap
        newstring = strrep(newstring, '___', '__'); 

        fprintf('    -> Removed tag: %s\n', tagToRemove);
    end
end

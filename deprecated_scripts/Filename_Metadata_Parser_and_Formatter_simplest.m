%% Filename Metadata Parser and Formatter
clear; clc;

% --- Configuration ---
inputDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD57\Standardized Testing\SHG Search\id33.2_50um_Sweep_tmp';  % Folder containing source files
outputDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD57\Standardized Testing\SHG Search\id33.2_50um_Sweep_new'; % Where the new files will go
inputExt  = '.txt';

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

files = dir(fullfile(inputDir, ['*', inputExt]));

for i = 1:length(files)
    oldName0 = files(i).name;
    newName = oldName0;
    newName = strrep(newName,'nm','');
    newName = strrep(newName,'P=','');
       
    
    % 4. File Operations
    fullInputPath  = fullfile(inputDir, oldName0);
    fullOutputPath = fullfile(outputDir, newName);
    
    % Using 'copyfile' instead of read/write since you mentioned "copy to a format"
    % This preserves the internal data exactly as is.
    [status, msg] = copyfile(fullInputPath, fullOutputPath);
    
    if status
        fprintf('Success: %s \n      -> %s\n', oldName0, newName);
    else
        fprintf('Error copying %s: %s\n', oldName0, msg);
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

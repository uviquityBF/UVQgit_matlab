%% Batch Internal Text Replacer (Recursive)
clear; clc;

% --- Configuration ---
rootFolder = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD57\Standardized Testing\SHG Search\'; 
targetExt  = '.txt';

% The strings you are looking for inside the files
oldString = 'slit=10um'; 
newString = 'slit=50um';

% SET TO true TO OVERWRITE FILES. Set to false to just "test" the search.
isWriteEnabled = true; 

% 1. Find all files recursively
% The '**/*.txt' pattern searches all subdirectories
filePattern = fullfile(rootFolder, ['**/*', targetExt]);
files = dir(filePattern);

fprintf('Found %d files. Searching for string: "%s"...\n', length(files), oldString);

matchCount = 0;

for i = 1:length(files)
    % Skip directories if they somehow get caught in the list
    if files(i).isdir, continue; end
    
    fullPath = fullfile(files(i).folder, files(i).name);
    
    % 2. Read the entire file content
    fileContent = fileread(fullPath);
    
    % 3. Check if the string exists
    if contains(fileContent, oldString)
        matchCount = matchCount + 1;
        fprintf('Match found in: %s\n', files(i).name);
        
        % 4. Perform the replacement
        updatedContent = strrep(fileContent, oldString, newString);
        
        % 5. Write back to file if enabled
        if isWriteEnabled
            fid = fopen(fullPath, 'w');
            if fid == -1
                error('Could not open file for writing: %s', fullPath);
            end
            fprintf(fid, '%s', updatedContent);
            fclose(fid);
        else
            fprintf('   [TEST MODE]: Replacement not saved.\n');
        end
    end
end

if isWriteEnabled
    fprintf('\nTask Complete. Updated %d files.\n', matchCount);
else
    fprintf('\nTask Complete. Found %d files containing the string (No changes saved).\n', matchCount);
end
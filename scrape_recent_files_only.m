%% File Extractor Script
% This script searches recursively for files matching a regex and 
% copies *THE MOST RECENT ONE of* them to a 'copies' directory at the top level.

function scrape_recent_files_only()
    % --- Configuration ---
%     rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\Run C';
%     rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\GaugeLots';
%    rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\';
%     rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\GaugeLots\GL02\';
    rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\Run G\AM164';
    % Pattern matches any file containing the string, regardless of extension
    
    pattern  = '_Overhead_Loss_Results - All'; 
    destDir  = fullfile(rootPath, 'copies_recent');

    if ~exist(destDir, 'dir'), mkdir(destDir); end

    % --- Step 1: Recursive search ---
    % We grab all files first to evaluate timestamps
    allFiles = dir(fullfile(rootPath, '**', '*.*'));
    allFiles = allFiles(~[allFiles.isdir]);

    % --- Step 2: Filter by Pattern ---
    matchIdx = ~cellfun(@isempty, regexp({allFiles.name}, pattern, 'once'));
    matchedFiles = allFiles(matchIdx);

    if isempty(matchedFiles)
        fprintf('No files matching the pattern were found.\n');
        return;
    end

    % --- Step 3: Group by Folder and Find Most Recent ---
    % We extract unique folders to find the newest file in each
    allFolders = {matchedFiles.folder};
    uniqueFolders = unique(allFolders);
    
    count = 0;
    for f = 1:length(uniqueFolders)
        currentFolder = uniqueFolders{f};
        
        % Skip the destination folder itself
        if contains(currentFolder, destDir), continue; end
        
        % Get all matched files in THIS specific folder
        thisFolderIdx = strcmp(allFolders, currentFolder);
        filesInThisFolder = matchedFiles(thisFolderIdx);
        
        % Find the file with the maximum (most recent) datenum
        [~, latestIdx] = max([filesInThisFolder.datenum]);
        latestFile = filesInThisFolder(latestIdx);
        
        % --- Step 4: Origin-Aware Copying ---
        sourcePath = fullfile(latestFile.folder, latestFile.name);
        [~, parentName] = fileparts(latestFile.folder);
        
        % Construct new name: ParentFolder_OriginalFileName
        newFileName = sprintf('%s_%s', parentName, latestFile.name);
        targetPath = fullfile(destDir, newFileName);

        [status, msg] = copyfile(sourcePath, targetPath);

        if status
            fprintf('Copied Latest: %s from %s (%s)\n', ...
                latestFile.name, parentName, latestFile.date);
            count = count + 1;
        else
            warning('Failed to copy %s: %s', latestFile.name, msg);
        end
    end

    fprintf('\nDone! Successfully copied %d most-recent files to %s\n', count, destDir);
end
%% File Extractor Script
% Recursively searches a folder tree for filenames matching a regex and
% copies matches into a destination folder at the top level. Two modes:
%   'all'                      - copies every match, renaming on collision
%                                 (formerly scrape_files_out_of_subfolders.m)
%   'most_recent_per_folder'   - groups matches by containing folder and keeps
%                                 only the newest (by datenum) file per folder,
%                                 prefixed with its parent folder name
%                                 (formerly scrape_recent_files_only.m)

function scrape_files()
    %% --- Configuration ---
    keep_mode = 'all';   % 'all' | 'most_recent_per_folder'
    rootPath  = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GLI\M3266\SD01\Standardized Testing\Overhead Loss';
    pattern   = '_Overhead_Loss_Results - All';   % regexp matched against filenames

    switch keep_mode
        case 'all'
            destDir = fullfile(rootPath, 'copies (all)');
        case 'most_recent_per_folder'
            destDir = fullfile(rootPath, 'copies_recent');
        otherwise
            error('Unknown keep_mode: %s', keep_mode);
    end

    %% --- Step 1: Create the destination folder ---
    if ~exist(destDir, 'dir')
        mkdir(destDir);
        fprintf('Created directory: %s\n', destDir);
    end

    %% --- Step 2: Search for all files recursively and filter by pattern ---
    allFiles = dir(fullfile(rootPath, '**', '*.*'));
    allFiles = allFiles(~[allFiles.isdir]);

    matchIdx = ~cellfun(@isempty, regexp({allFiles.name}, pattern, 'once'));
    matchedFiles = allFiles(matchIdx);

    if isempty(matchedFiles)
        fprintf('No files matching the pattern were found.\n');
        return;
    end

    %% --- Step 3: Copy per mode ---
    switch keep_mode
        case 'all'
            count = copy_all(matchedFiles, destDir);
        case 'most_recent_per_folder'
            count = copy_most_recent_per_folder(matchedFiles, destDir);
    end

    fprintf('\nDone! Successfully copied %d file(s) to %s\n', count, destDir);
end


%% ================= Helper functions =================

function count = copy_all(matchedFiles, destDir)
% Copies every matched file into destDir, renaming on filename collision.
    count = 0;
    for i = 1:numel(matchedFiles)
        % Skip files already inside the destination folder to avoid infinite loops
        if contains(matchedFiles(i).folder, destDir)
            continue;
        end

        fileName = matchedFiles(i).name;
        sourcePath = fullfile(matchedFiles(i).folder, fileName);
        targetPath = fullfile(destDir, fileName);

        if exist(targetPath, 'file')
            [~, name, ext] = fileparts(fileName);
            targetPath = fullfile(destDir, sprintf('%s_%d%s', name, count, ext));
        end

        [status, msg] = copyfile(sourcePath, targetPath);
        if status
            fprintf('Copied: %s\n', fileName);
            count = count + 1;
        else
            warning('Failed to copy %s: %s', fileName, msg);
        end
    end
end


function count = copy_most_recent_per_folder(matchedFiles, destDir)
% Groups matched files by containing folder and copies only the newest
% (by datenum) file per folder, prefixed with its parent folder name.
    count = 0;
    allFolders = {matchedFiles.folder};
    uniqueFolders = unique(allFolders);

    for f = 1:numel(uniqueFolders)
        currentFolder = uniqueFolders{f};

        % Skip the destination folder itself
        if contains(currentFolder, destDir)
            continue;
        end

        filesInThisFolder = matchedFiles(strcmp(allFolders, currentFolder));
        [~, latestIdx] = max([filesInThisFolder.datenum]);
        latestFile = filesInThisFolder(latestIdx);

        sourcePath = fullfile(latestFile.folder, latestFile.name);
        [~, parentName] = fileparts(latestFile.folder);
        newFileName = sprintf('%s_%s', parentName, latestFile.name);
        targetPath = fullfile(destDir, newFileName);

        [status, msg] = copyfile(sourcePath, targetPath);
        if status
            fprintf('Copied Latest: %s from %s (%s)\n', latestFile.name, parentName, latestFile.date);
            count = count + 1;
        else
            warning('Failed to copy %s: %s', latestFile.name, msg);
        end
    end
end

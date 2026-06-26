function overheadLoss_open_CSVs_and_compile3()
    clear
    close all;

    outFile = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\TMP.csv';

    % Remove old output file and initialize
    if exist(outFile, "file")
        delete(outFile);
    end

    % Directory to scan
%    rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run A\2025_04_29_RunA_OverheadLoss';
%     rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Round2 (xiyuan)\2025_04_09_Round2_Loss\Overhead Loss';
%     rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run A\RunA Survey';
%     rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136\SD43\R&R Study#3\Overhead Loss';
%     rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136\SD43\R&R Study Follow Up\Overhead Loss' 
%    rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136\SD43\R&R Study #1\Overhead Loss';
%    rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run A\Repeatability Study (RunA_RunB)\Results\Perry';
%     rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run A\Repeatability Study (RunA_RunB)\Results\Jacob';
%    rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run A\Repeatability Study (RunA_RunB)\Results\Brent';

  rootDir = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run C\AM123\'
    
    
    % Find all CSV files
    fileStruct = dir(fullfile(rootDir, '**', '*Params*.csv'));
    N = length(fileStruct);

    fprintf("Found %d CSV files.\n", N);

    % Storage
    allHeaders = string([]);     % master header list
    allRows    = {};             % each row will be a containers.Map

    % ---- Read each CSV --------------------------------------------------
    for k = 1:N
        filePath = fullfile(fileStruct(k).folder, fileStruct(k).name);
        fprintf("Processing: %s\n", fileStruct(k).name);

        [headers, values] = readCSV_headers_and_values(filePath);
        s = split(filePath,'\');        filename = s{end};
        ss = split(filePath,s{end});    path = ss{1};            
        values(end+1) = path;
        values(end+1) = filename;
        headers(end+1) = "csv_path";
        headers(end+1) = "csv_filename";
        
        
        % store dynamic row in a map
        rowMap = containers.Map(headers, values);

        % update master header list
        allHeaders = union(allHeaders, headers, "stable");

        allRows{end+1} = rowMap;
    end

    % ---- Write unified table to TMP.csv --------------------------------
    writeMasterCSV(outFile, allHeaders, allRows);

    fprintf("DONE.  Wrote output to:\n%s\n", outFile);
end




%% ========================================================================
function [headers, values] = readCSV_headers_and_values(filepath)
    % Read CSV using detectImportOptions, convert to strings
    opts = detectImportOptions(filepath);
    T = readtable(filepath, opts);
    C = string(table2cell(T));

    headers = C(:,2);    % second column = header names
    values  = C(:,3);    % third column  = data to place under those headers
end



%% ========================================================================
function writeMasterCSV(outFile, allHeaders, allRows)
    % Open file
    fID = fopen(outFile, 'w');

    % ---- Write header row ----
    headerLine = strjoin(allHeaders, ',');
    fprintf(fID, "%s\n", headerLine);

    % ---- Write each data row ----
    for r = 1:length(allRows)
        rowMap = allRows{r};

        % Preallocate row with blanks
        rowVals = strings(1, length(allHeaders));

        % Fill in values where they exist
        for h = 1:length(allHeaders)
            key = allHeaders(h);
            if isKey(rowMap, key)
                rowVals(h) = rowMap(key);
            end
        end

        % Write row
        fprintf(fID, "%s\n", strjoin(rowVals, ','));
    end

    fclose(fID);
end
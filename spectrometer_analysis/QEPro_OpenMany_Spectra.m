% Open QEPro Data Files
% export single matrix CSV
%
% Set file_list_csv to a CSV path (one filepath per line) to load the file
% list from that CSV instead of the interactive uigetfile dialog.
function QEPro_OpenMany_Spectra()
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\spectra');
%% get data
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run C\SD33 (TD5)\SHG Search\id11.2\WL Sweep'
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Tools\Slit Width (QE Pro)\2025_06_09__TiSaph ND20+03 + FiberG';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\2025_06_30_gc310_Pulsed';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run C\SD41 (TD4)\SHG Search\id31.2\2025_06_05 Fine WL Sweep - 225nm alignment - Z optimized';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_28_LabData\CAL_10umFiber\spectra v Z';

    % Set to a CSV path (one filepath per line) to bypass uigetfile and load
    % the file list from there instead.
    file_list_csv = '';

%% get list of files
    if isempty(file_list_csv)
        cd(start_path)
        [filenames,path] = uigetfile('*.txt','Multiselect','on');
        Nf = length(filenames);
        if isstr(filenames)
            filename = filenames; clear filenames;
            filenames{1} = filename;
            Nf=1;
        end
        filepaths = cell(Nf,1);
        for kf=1:Nf
            filepaths{kf} = [path,filenames{kf}];
        end
        out_dir = path;
    else
        filepaths = read_filepath_list(file_list_csv);
        Nf = length(filepaths);
        [out_dir_raw,~,~] = fileparts(file_list_csv);
        out_dir = [out_dir_raw,'\'];
    end

% Loop Over Files
    for kf=1:Nf
        clear DAT;
        try
            DAT = read_spectrum_file(filepaths{kf}, 'oceanoptics');
        catch
            disp(filepaths{kf}); % couldn't open/parse -- skip
            continue;
        end
        M(:,kf) = DAT.counts; %copy spectrum to matrix
        X(:,kf) = DAT.wl;
        Tint(kf,1) = DAT.Tint;
    end

    csvwrite([out_dir,'allspectra.csv'],[X(:,1),M]);
    disp(filepaths');

end

function C = read_filepath_list(filepath)
    fID = fopen(filepath);
    k=1;
    C = {};
    while ~feof(fID)
        C{k,1} = fgetl(fID); k=k+1;
    end
    fclose(fID);
end

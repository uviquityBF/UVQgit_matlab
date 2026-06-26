% Open QEPro Data Files
% export single matrix CSV
%
function QEPro_OpenMany_Spectra()
%% get data
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    
%get list of files from CSV
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples';
    filename = 'tmp.csv';
    files_list_csv = [dirs.d,'\',filename];
    
    filepath_and_names = read_filepath_list(files_list_csv );
    Nf = length(filepath_and_names);
    
% Loop Over Files
    for kf=1:Nf
        clear DAT;
        DAT  =  readQEPro_SpectrumFile([filepath_and_names{kf}]);
        if isstruct(DAT)
            M(:,kf) = DAT.A(:,2); %copy spectrum to matrix
            X(:,kf) = DAT.A(:,1);
            Tint(kf,1) = DAT.Tint;
        if Nf <= 20
            figure; plot(DAT.A(:,1),DAT.A(:,2)); hold all; plot(DAT.A(isel,1),DAT.A(isel,2));
            title({replace(filenames{kf},'_',' '); ['sig = ',num2str(signal(kf,1))]});
        end
        end
        
    end

    csvwrite([dirs.d,'\allspectra.csv'],[X,M]);
    disp(filenames');

end

function C = read_filepath_list(filepath)
    fID = fopen(filepath);
    k=1;
    while ~feof(fID)
        C{k,1} = fgetl(fID); k=k+1;
    end
    fclose(fID)
    
%     % Read CSV using detectImportOptions, convert to strings
%     opts = detectImportOptions(filepath);
%     T = readtable(filepath, opts);
%     C = string(table2cell(T));
end




 %% Load QE Pro File

function DAT  =  readQEPro_SpectrumFile(filepath_and_name)
    fID = fopen(filepath_and_name);
    clear DAT;
    if fID~=-1

        %Skip Header
        Nskip = 14;
        frewind(fID);     	
        for kk=1:Nskip
            tline1 = fgetl(fID);
            [newstrings,matches] = split(tline1,'(sec):');
            if length(matches)>0
                DAT.Tint = str2num(newstrings{2});
            end
        end

        %read spectral data
         kk=1;
        while feof(fID)==0  %while NOT end of file
            tline1 = fgetl(fID);

            %[1] split values by delimiter
            %               [vals,matches] = split(tline1,'\t');  %DOESNT WORK
            %[2] split values by 8 characters for each %DOESNT WORK
            %             val1 = num2str(tline1(1:7));
            %             val2 = num2str(tline1(9:end));
            %[3] split using "sscanf"
            vals = sscanf(tline1,'%f') ; %works
            if length(vals)>=2
                wl =    vals(1);
                counts =  vals(2);
                DAT.A(kk,1) =   wl;
                DAT.A(kk,2) =  counts;
                kk=kk+1;
            end
        end

%         figure; plot(DAT(:,1),DAT(:,2));
    else
        DAT=0;
        disp(filepath_and_name)
    end %if fid~=1
    
end

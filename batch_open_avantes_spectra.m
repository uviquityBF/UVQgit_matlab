%% Open Avantes Spectra (ASCII)
function batch_open_avantes_spectra()
    clear;
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_04_21_LabData\PulseDiode 2025_04_24\spectra';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_04_21_LabData\PulseDiode 2025_04_23\spectra';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_05_05_LabData';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_05_05_LabData\ECL (cw) -- gc309,310';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_05_05_LabData\pulsed diode - dfb383 (LIV)';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\DiodeLasers\2025_06_30_gc310_Pulsed';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_06_30_LabData\Q - wideband spectrum repeatability';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_14_LabData\Q - using Spectrometer\spectra\Round0_ALN1C';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_14_LabData\Q - using Spectrometer\spectra\Round2_SGH3B1s';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_08_04_LabData\Q_v_Temperature\RLB#1';
    
    %Select Files
    cd(dirs.d);
    [filenames,start_path] = uigetfile({'*.txt'},'Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end

    %Open Each File
    for kf=1:Nf
        disp(filenames{kf});
        A = dlmread([start_path,'/',filenames{kf}],';',9,0);
        DAT.I(:,kf)     = A(:,2);
        DAT.Dark(:,kf)  = A(:,3);
        DAT.Ref(:,kf)   = A(:,4);
        if kf==1;
            WL              = A(:,1);
        end
        labels{kf} = replace(filenames{kf},'_',' ');
    end
    

    %output
    out = [WL,DAT.I];
    csvwrite('temp.txt',out);
    
    disp(filenames');

        %show data
%     DAT.I(:,end-1:end) = DAT.I(:,end-1:end)/10;
    figure; plot(WL,DAT.I); legend(labels);

    
    
    
end



 


  




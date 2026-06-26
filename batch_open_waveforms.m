%% Open Waveforms (ASCII)
function batch_open_waveforms()
    clear;
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_04_28_LabData\2025_04_29_diode_pulse_l450p1600mm__Waveform-v-Grating\waveforms_equal_grating_steps';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_05_12_LabData\diode pulse_ECL gc304\vary_WL_with_OscopeData\waveforms(oscope)';
    %Select Files
    cd(dirs.d);
    [filenames,start_path] = uigetfile({'*.csv'},'Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end

    %Open Each File
    for kf=1:Nf
        disp(filenames{kf});
        A = csvread([start_path,'/',filenames{kf}],44,0);
        DAT.I(:,kf)     = A(:,2);
        if kf==1;
            t           = A(:,1);
        end
        labels{kf} = replace(filenames{kf},'_',' ');
    end
    

    %output
    out = [t,DAT.I];
    csvwrite('temp.txt',out);
    
    disp(filenames');

        %show data
    figure; plot(t,DAT.I); legend(labels);

    
%     figure; plot(t,DAT.I(:,1),...
%             t,DAT.I(:,2),...
%             t,DAT.I(:,3),...
%             t-2.4e-8,DAT.I(:,4),...
%             t,DAT.I(:,5),...
%             t,DAT.I(:,6)); 
%     legend(labels);

    
    
end



 


  




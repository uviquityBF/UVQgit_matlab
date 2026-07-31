% Open QEPro Data Files
%
%
function MAIN()
%% get data
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\spectra');
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    filepath = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    filepath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD57\Standardized Testing\SHG Search\id1.2_50um_Sweep'
%     filepath = [uigetdir(start_path),'/'];

    cd(start_path)
    [filenames,path] = uigetfile('*.txt','Multiselect','on');
     Nf = length(filenames);
     if isstr(filenames)
         filename = filenames; clear filenames;
         filenames{1} = filename;
         Nf=1;
     end


% Loop Over Files
    for kf=1:Nf

        DAT  =  read_spectrum_file([path,filenames{kf}], 'oceanoptics');
        M(:,kf) = DAT.counts; %copy spectrum to matrix
        BG(kf,1) = mean( DAT.counts(end-10:end) );
        BGsigma = std(DAT.counts(end-10:end) );
        BGrange = [DAT.wl(end-10),DAT.wl(end)];

        thresh = 200;   %BG(kf,1) + 50*BGsigma
        isel = find(DAT.counts > thresh);
        signal(kf,1) = sum(DAT.counts(isel));
        numbins(kf,1) = length(isel);
        Tint(kf,1) = DAT.Tint;
        if Nf <= 20
            figure; plot(DAT.wl,DAT.counts); hold all; plot(DAT.wl(isel),DAT.counts(isel));
            title({replace(filenames{kf},'_',' '); ['sig = ',num2str(signal(kf,1))]});
        end
    end
    tmp = split(filenames{kf},'QEP064681');
    titlestring = replace(tmp{1},'_',' ');
    figure; imagesc(M);  colorbar;   title({titlestring;'linear scale'});
    colorbar;    ylabel('wavelength[nm]');  xlabel('spectrum count');
    figure;  imagesc([1:size(M,2)],DAT.wl,abs(log10(M))); title({titlestring;'log10 scale'});
    colorbar;    ylabel('wavelength[nm]');  xlabel('spectrum count');

    Mavg  = mean(M,2);
    figure; plot(DAT.wl,Mavg);  xlabel('wavelength [nm]');  title({titlestring;'linear scale'});
    grid on;  ylim([50,150]*15); xlim([250,400]);



    csvwrite([path,'temp.csv'],[signal,numbins,BG,Tint]);
    csvwrite([path,'allspectra.csv'],[DAT.wl,M]);
    disp(filenames');

end

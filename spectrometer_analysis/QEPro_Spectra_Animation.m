% 1. Open Spectra Files 
 
%
function QEPro_Spectra_Animation()
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\spectra');
    close all;
    clear;
    program_name = 'image analysis';
    
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_23_LabData\';
    
%% select SPECTRA STACK
   
    %Select by Hand
    cd(start_path);
    [filenames,path] = uigetfile('*.txt','Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end
%% Read All SPECTRA Files
    for kf = 1:Nf
       DAT{kf} =  read_spectrum_file([path,filenames{kf}], 'oceanoptics')
    end
    
%% select BG SPECTRUM
    %Select by Hand
    cd(start_path);
    [BG_filename,BG_path] = uigetfile('*.txt','Multiselect','off');
    if BG_filename==0   BGflag = 0;
    else                BFflag = 1;
    end
    
    if isstr(BG_filename)   BG_filename = BG_filename;  
    end
%% Read BG SPECTRUM
    BG =  read_spectrum_file([BG_path,BG_filename], 'oceanoptics');
    
    
s = split(filenames{1},'set') 
s2 = split(s{2},'__') 
WL = str2num(s2{1});

%% Generate Animation
%set up video
v = VideoWriter('test3.avi');
v.FrameRate = 10;
open(v);
hf = figure;
for kf = 1:Nf
    s = split(filenames{kf},'set');
    s2 = split(s{2},'__') ;
    pumpWL = str2num(s2{1});

    %figure update
    figure(hf);
    WL = DAT{kf}.wl;
    Y = DAT{kf}.counts - BG.counts;
    plot(WL,Y); xlim([280,340]);
    title(['pump WL=',num2str(pumpWL),'nm   frame=',num2str(kf)]);
    ylim([0,40]);
    
    %generate video
    movie_F(kf) = getframe(gcf);
    writeVideo(v,movie_F(kf));
    pause(0.25);    
end

close(v);

movie(gcf,movie_F,10,1);


end
function signal_inside_radius()
    clear;
    close all;
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    
    
    
    %........ DIRECTORIES where relevant data is  
    dirs.input = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_07_LabData\straylight';
    DoTwoCircles = 1;
   
    
     dirs.input = [dirs.input, ''];
     cd(dirs.input);
      
    %select files to plot
    [IMGfile_csv,IMGpath,dummy] = uigetfile({'*.csv'},'Select DATA File to ','MultiSelect','off');
    
    %read data and find peak
    A = csvread([IMGpath,IMGfile_csv]);
    [maxval,imax] = max(A(:)); 
    [ymax,xmax]=ind2sub(size(A),imax);    
    %define R from main peak
    [XX,YY]=meshgrid([1:size(A,2)],[1:size(A,1)]);
    RR = ( (XX-xmax).^2 + (YY-ymax).^2 ).^0.5;    
    %define R from secondary point
    if DoTwoCircles == 1
        Yoffset = 175; %pixels
        RR2 = ( (XX-xmax).^2 + (YY-ymax-Yoffset).^2 ).^0.5;
    else RR2 = RR; end
    figure; imagesc(RR)
    figure; imagesc(RR2)
    
    %Vary Radius
    R_limit_min = 5;
    R_limit_max = 54;
    R_limit_list = linspace(R_limit_min,R_limit_max,50)';
    for k=1:length(R_limit_list)
        R_limit = R_limit_list(k);
        isel1 = find(RR(:)<R_limit);
        if DoTwoCircles == 1 
                isel2 = find(RR2(:)<R_limit);
        else isel2 = isel1; end
        isel = unique([isel1;isel2]); 
        signal_v_R(k,1) = sum( A(isel) );        
    end;
    figure; plot(R_limit_list,signal_v_R/signal_v_R(1)); grid on;
    ylabel(['signal(R_l_i_m_i_t) / signal(R_l_i_m_i_t=',num2str(R_limit_min),')']);
    xlabel('R_l_i_m_i_t [pixels]');
    title(replace(IMGfile_csv,'_',' '));
    
    disp(signal_v_R/signal_v_R(1));
    
    %Illustration: analyze single Radius
    R_limit = R_limit_max;
        isel1 = find(RR(:)<R_limit);
        if DoTwoCircles == 1 
                isel2 = find(RR2(:)<R_limit);
        else isel2 = isel1; end
        isel = unique([isel1;isel2]); 
    signal_sel = sum( A(isel) );
    
    thetas=linspace(0,2*pi,100)';
    hf = figure;
    imagesc(log10(A));
    hold all; plot(R_limit*cos(thetas)+xmax,R_limit*sin(thetas)+ymax,'Color','k','Linewidth',2);
    title(['Signal inside circle: ' num2str(round(signal_sel,0))]);
    if DoTwoCircles == 1 
     hold all; plot(R_limit*cos(thetas)+xmax,R_limit*sin(thetas)+ymax+Yoffset,'Color','k','Linewidth',2);    
    end
    
    
    
    
    

function DLPro_ModeMap()
%% get data
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\Konrad_Testing\Data\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_07_LabData\modemap repeatability';
%     filepath = [uigetdir(start_path),'/'];

    %get file names
   cd(start_path)
   [filename1,path] = uigetfile('*.csv','Multiselect','off')    
   A = csvread([path,filename1],1,0);
   [mm1.VV,mm1.II, mm1.LAM, mm1.I,mm1.V] = make_modemap(A,[],[])
   title(replace(filename1,'_',' '))

   [filename2,path] = uigetfile('*.csv','Multiselect','off')    
   B = csvread([path,filename2],1,0);
   [mm2.VV,mm2.II, mm2.LAM, mm2.I ,mm2.V] = make_modemap(B,mm1.I,mm1.V);
   title(replace(filename2,'_',' '))
   
   deltaLAM = (mm1.LAM - mm2.LAM); 
   figure; imagesc(mm1.V, mm1.I, abs(deltaLAM)); colorbar;  caxis([0,0.1]);
   title(['difference: ',replace(filename1,'_',' '),' - ', replace(filename2,'_',' ')]);
   
   %COLOR SLICES 
   LAM = mm1.LAM; V = mm1.V; I = mm1.I;
   wldel = 0.15; %[nm]
   list.wlctr = [449 : wldel : 451]';
   list.wlctr = [450.15 , 450.39]';
   for k=1:length(list.wlctr) 

       wlctr =  list.wlctr(k);  %[nm]

       % apply filter
       isel = find( wlctr-wldel/2 <= LAM(:) & LAM(:) <=  wlctr+wldel/2 );
       %[rsel,csel] = ind2sub(size(LAM),isel);
       Nsel = length(isel)
       LAMsel = NaN*ones(size(LAM));
       LAMsel(isel) = LAM(isel);

       %plot
       figure;  
       imagesc(V,I,flipud(LAMsel));
       colorbar(); 
       title(['WL_C_T_R = ',num2str(wlctr),'nm']);
       
   end
   
%    % image segmentation
%    imageSegmenter(LAM)
%    
%    level = multithresh(LAM,20)
   
end
 
function [VV,II,LAM, I,V] = make_modemap(A,I,V)
 %  Current_mA	Piezo_V     Frequency_THz	Power_uW	Temperature_C	Pressure_mbar
   raw.I = A(:,1);
   raw.V = A(:,2);
   raw.f = A(:,3);
   raw.T = A(:,4);
   raw.prs = A(:,5);
    
    
   %histogram of all WLz 
   raw.WL           = 3e8 ./(raw.f*1e12)  *1e9;
   binsize_nm       = 0.001;
   xbins            = [449:binsize_nm:452]';
   [hvals,xvals]    = hist(raw.WL,xbins);
   figure;          plot(xvals,hvals);
   
   % sort through data
   if length(I)==0    I = unique(raw.I);  end
   if length(V)==0    V = unique(raw.V);  end
   
   [VV,II] = meshgrid(V, I);
   LAM = griddata(raw.V,raw.I,raw.WL,VV,II);   %[nm]
   figure; imagesc(V,I,LAM);  colorbar();  caxis([448,452]);
%    LAM = 3e8 ./(griddata(raw.V,raw.I,raw.f,VV,II)*1e12)  *1e9 ;   %[nm]
   
end

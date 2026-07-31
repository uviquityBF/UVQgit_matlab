
%% Simulate combination of a narrow and wide slit with narrow and wide spectrum
%  Method:  Convolve narrow and wide guassian with narrow and wide boxcar 


%%

clear; 
Xmax = 0.3;           %[mm]
dX   = 0.0001;      %[mm]
x01 = 0;        
x02 = 0;
width1_um = 10;     %[um]   10um Slit
width2_um = 200;    %[um]   10um Slit
fwhm1_um = 1.5;     %[um]   0.1nm linewidth: equivalent slit width
fwhm2_um = 30;      %[um]   2nm linewidth: equivalent slit width


%% define functions
close all;
width1 = width1_um/1000;
width2 = width2_um/1000;
fwhm1 = fwhm1_um / 1000;
fwhm2 = fwhm2_um / 1000;

x = [-Xmax:dX:Xmax]';   %units = [mm]

yslit1 = double( abs(x-x01) < 0.5*width1 );
yslit2 = double( abs(x-x02) < 0.5*width2 );

yg1 = exp(-(x-x01).^2 / (2*fwhm1)^2 );
yg2 = exp(-(x-x02).^2 / (2*fwhm2)^2 );


% figure; plot(x,[yslit1,yslit2]);
% figure; plot(x,[yg1,yg2]);
figure; 
subplot(2,1,1); plot(x,[yslit1,yslit2]);title('Input Assumptions');
legend(['slit1: ',num2str(width1_um),'um'],['slit2: ',num2str(width2_um),'um']);
subplot(2,1,2); plot(x,[yg1,yg2]);  xlabel('position on spectrometer detector array');  
legend('spectral line1 (yg1)','spectral line2 (yg2)');


%% convolve
yg1_slit1 = conv(yg1,yslit1,'same');
yg2_slit1 = conv(yg2,yslit1,'same');
yg1_slit2 = conv(yg1,yslit2,'same');
yg2_slit2 = conv(yg2,yslit2,'same');

figure; plot(x,[yg1_slit1,yg1_slit2]);
legend(['yg1 slit1: ',num2str(width1_um),'um'],...
       ['yg1 slit2: ',num2str(width2_um),'um']); 
xlabel('position on spectrometer detector array');    
title(['Q Tune comparison:  linewidth = 0.1nm']);


figure; plot(x,[yg2_slit1,yg2_slit2]);
legend(['yg2 slit1: ',num2str(width1_um),'um'],...
       ['yg2 slit2: ',num2str(width2_um),'um']); 
title(['TiSapphire comparison:  linewidth = 2nm']);
xlabel('position on spectrometer detector array');    



% figure; plot(x,[yg1_slit1,yg2_slit1,yg1_slit2,yg2_slit2]);
% legend('yg1 slit1','yg2 slit1', 'yg1 slit2', 'yg2 slit2');




















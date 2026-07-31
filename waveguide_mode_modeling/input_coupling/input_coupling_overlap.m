% Input Coupling Overlap Calculation

%[1]  get PSF
%[2]  define ellipsoid of mode
%[3] calculate overlap
function input_coupling_overlap()
    close all;
    clear;

    dirs.d1 = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024\2024_05_01-02\2024_05_02_InputSpot1_for_HDR';
    dirs.d2 = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024\2024_05_01-02\2024_05_02_InputSpot2_for_HDR';
    dirs.d = dirs.d1;

    %% get PSF
    scale_um_per_pixel = (40 / 360 * (50/100));
    
    cd(dirs.d);
    [filename,path] = uigetfile();
    A = load([path,filename]);
    [max_val,imax]  = max(A.I_HDR(:));
    [Nr,Nc]=size(A.I_HDR);
    [max_r,max_c] = ind2sub([Nr,Nc],imax);
    Xpix = [1:Nc] - max_c;
    Ypix = [1:Nr]' - max_r;
    figure; imagesc(Xpix,Ypix,A.I_HDR); axis equal; xlim([-100,100]); ylim([-40,40]);
    [XXpix,YYpix] = meshgrid(Xpix,Ypix); 
  %  figure; imagesc(Xpix,Ypix,XXpix);
   % figure; imagesc(Xpix,Ypix,YYpix);

    %% define Gaussien
    halfwidth_um.x = 1.25;      %units = pixels
    halfwidth_um.y = 0.175;       %units = pixels

    halfwidth_pix.x = halfwidth_um.x/scale_um_per_pixel;      %units = pixels
    halfwidth_pix.y = halfwidth_um.y/scale_um_per_pixel;       %units = pixels
    Z = (1/sqrt(2*halfwidth_pix.x)).*exp( -0.5*(XXpix/(halfwidth_pix.x)).^2 ) ...
        .* (1/sqrt(2*halfwidth_pix.y)).*exp( -0.5*(YYpix/(halfwidth_pix.y)).^2 );    
    figure; imagesc(Xpix,Ypix,Z); axis equal; xlim([-100,100]); ylim([-40,40]);


    %% calculate overlap
    overlap_matrix  = ( A.I_HDR.*Z ) / sum(A.I_HDR(:))*sum(Z(:));
    overlap = sum(overlap_matrix(:));

end


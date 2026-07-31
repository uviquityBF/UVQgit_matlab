% Index Surfaces (Analytic)

function  null = modal_index_surfaces()
    clear;
    close all;
    WL      = 600;
    wdt     = 700;
       
%     %FUND dependence
%     a1 = -0.001563;     %WL dependence
%     b1 =  0.00070;      %width dependence
%     f1 =  0;            %height dependence
%     a2 = 1.07E-06;      %WL dependence ^2
%     b2 = -4.92E-07;     %width dependence ^2
%     f2 = 0;
% 
%     %SHG Dependence
%     c1 = -0.01054;
%     d1 = 0.0027;  
%     g1 =  0;
%     c2 = 1.21E-05;
%     d2 = -1.40E-06;
%     g2 = 0.00E+00;
     
   
    % Fundamental -- Fit to Width 
    n_f0.w      = 1.93893;     %1.69129401444789;     %fit to width
    WL_f0.w     = 300;
    wdt_f0.w    = 700;
    b1          = 9.57e-6;      %6.98e-4    %width dependence
    b2          = -4.92E-07;    %-4.92E-07; %width dependence ^2
    
    % Fundamental -- Fit to WL
    n_f0.wl     = 1.941420378;  %2.3133;    %fit to WL
    WL_f0.wl    = 300;
    wdt_f0.wl   = 700;
    h_f0.wl     = 390;
    a1          = -0.0009252864146;     %-0.001563;     %WL dependence
    a2          = 1.06E-06;            %1.07E-06;      %WL dependence ^2
    	 
    % SHG (TM40) -- Fit to Width    
    n_s0.w      = 1.94;         %0.751;     %fit to width
    WL_s0.w     = 300;
    wdt_s0.w    = 700;
    h_s0.w      = 390;
    d1          = 7.26e-4;      %0.00268; 
    d2          = -1.4e-6;      %-1.40E-06;    
    
    % SHG (TM40) -- Fit to WL 
    n_s0.wl      = 1.94107;      %3.979233;    
    WL_s0.wl     = 300;
    wdt_s0.wl    = 700;
    h_s0.wl      = 390;
    c1          = -0.003427;   %-0.01054;    
    c2          = 1.33e-5;      %1.21E-05;    

%     % SHG (TM02) -- Fit to Width    
%     n_s0.w      = 1.930527554;      %fit to width
%     WL_s0.w     = 300;
%     wdt_s0.w    = 700;
%     h_s0.w      = 390;
%     d1          = -3.75E-05;       
%     d2          = -3.17E-07;    
%     
%     
%     % SHG (TM02) -- Fit to WL 
%     n_s0.wl      = 1.923796429;        
%     WL_s0.wl     = 300;
%     wdt_s0.wl    = 700;
%     h_s0.wl      = 390;
%     c1          = -0.003345714286;   
%     c2          = 1.52E-05;      
   

    %plot surfaces
    vec_WL = [200:10:500]';
    vec_wdt = [300:10:700]';
    
    vec_n_f.wl = n_f0.wl + a1*(vec_WL-WL_f0.wl) + a2*(vec_WL-WL_f0.wl).^2;
    vec_n_s.wl = n_s0.wl + c1*(vec_WL-WL_s0.wl) + c2*(vec_WL-WL_s0.wl).^2;
    figure; plot(vec_WL, [vec_n_f.wl,vec_n_s.wl]); xlabel('wavelength'); title('n v WL @ width=700')
        
    vec_n_f.w = n_f0.w + b1*(vec_wdt-wdt_f0.w) + b2*(vec_wdt-wdt_f0.w).^2;
    vec_n_s.w = n_s0.w + d1*(vec_wdt-wdt_s0.w) + d2*(vec_wdt-wdt_s0.w).^2;
    figure; plot(vec_wdt, [vec_n_f.w,vec_n_s.w]); xlabel('width'); title('n v width @ WL=300')
    
    
    %generate surfaces    
    [XX,YY] = meshgrid(vec_WL,vec_wdt);
    
    NN_f.w = n_f0.w + a2*(XX-WL_f0.w).^2 + a1*(XX-WL_f0.w) + b2*(YY-wdt_f0.w).^2 + b1*(YY-wdt_f0.w);
    NN_f.wl = n_f0.wl + a2*(XX-WL_f0.wl).^2 + a1*(XX-WL_f0.wl) + b2*(YY-wdt_f0.wl).^2 + b1*(YY-wdt_f0.wl);
    
    NN_s.w = n_s0.w + c2*(XX-WL_s0.w).^2 + c1*(XX-WL_s0.w) + d2*(YY-wdt_s0.w).^2 + d1*(YY-wdt_f0.w);
    NN_s.wl = n_s0.wl + c2*(XX-WL_s0.wl).^2 + c1*(XX-WL_s0.wl) + d2*(YY-wdt_f0.wl).^2 + d1*(YY-wdt_f0.wl);

    %plot - fundamental
    figure; plot3([vec_WL],[wdt_f0.wl*ones(length(vec_WL),1)],[vec_n_f.wl],'Linewidth',4); grid on;
    hold all; plot3([WL_f0.w *ones(length(vec_wdt),1)],[vec_wdt],[vec_n_f.w],'Linewidth',4); grid on;
    hold on; mesh(XX,YY,NN_f.w); hold all; mesh(XX,YY,NN_f.wl); title('fund');

    %plot - shg
    figure; plot3([vec_WL],[wdt_s0.wl*ones(length(vec_WL),1)],[vec_n_s.wl],'Linewidth',4); grid on;
    hold all; plot3([WL_s0.w *ones(length(vec_wdt),1)],[vec_wdt],[vec_n_s.w],'Linewidth',4); grid on;
    hold on; mesh(XX,YY,NN_s.w); hold all; mesh(XX,YY,NN_s.wl); title('SHG');

    
    
    
    figure; surf(XX,YY,NN_f.w); hold all; surf(XX,YY,NN_s.w);
    xlabel('WL'); ylabel('width');

    figure; surf(XX,YY,NN_f.wl); hold all; surf(XX,YY,NN_s.wl);
    xlabel('WL'); ylabel('width');
    
    
    figure; surf(XX,YY,abs(NN_f.w-NN_s.w)); 
        xlabel('WL'); ylabel('width');
    

    for k = 1:length(vec_wdt)
        wdt = vec_wdt(k);
        A = (a2-c2);
        B = -2*(a2*WL_f0.w - c2*WL_s0.w) + (a1 - c1);
        C = a2*WL_f0.w^2 - c2*WL_s0.w^2 + c1*WL_s0.w - a1*WL_f0.w;
        K = ( d2*(wdt-wdt_s0.w)^2 + d1*(wdt-wdt_s0.w)  ) - ...
                ( b2*(wdt-wdt_f0.w)^2 + b1*(wdt-wdt_f0.w)  );
        
        WL_root(k,1) = -B/(2*A) + (1/(2*A))*sqrt(B^2 - 4*A*(C-K) );
    end
    figure; plot(WL_root,vec_wdt); xlim([200,500]); ylim([300,700]);
    
end




















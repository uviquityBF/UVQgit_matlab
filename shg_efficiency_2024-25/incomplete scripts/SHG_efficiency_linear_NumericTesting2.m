

function SHG_efficiency_linear_NumericTesting2()
    c.c0                            = 3e8; %[m/s]
    c.eps0                          = 8.85E-12;  %[C/(V.m)]
    c.dz_um                         = 0.25;     %[um]  stepsize for numeric integration 

    %parameters define:
    modn.Pp                     = 0.1;      %[W]
    modn.L_waveguide_um         = 2000;   %[um]
    modn.Loss_pump_dBcm         = 0;
    modn.Loss_SH_dBcm           = 100;      
    modn.Abs_SH_dBcm            = 0;
    modn.lam1                   = 450;      %[nm]
    modn.d_xx                   = 4.7;      %[pm/V]
    modn.n1                     = 2.15;     %[um]
    modn.n2                     = 2.6;      %[um]
    modn.x_waveguide_um         = 0.3365;   %[um]
    modn.y_waveguide_um         = 0.381;    %[um]
    modn.overlap_frac           = 0.03;
    modn.dispersion_slope_diff  = 0.0089;   %[1/nm]
    modn.modeindex_crossing     = 1.94;     %unitless
    
    % SINGLE RUN
    modn.delta_lam              = 0.1;    %this is delta from ideal PhaseMatch WL.  it is sum of:    
%     modn.Loss_SH_dBcm           = 100;   
%     modn.Abs_SH_dBcm            = 100; 
    [tmpmodn,PscattRatio,P3_guided_end,P3b_scatt_end, P3scatt_max] = integrate_z(c,modn,1);    
    
    
    
    %SWEEP: DETUNE & ALPHA3
    sweepvals_delta_lam     = [0.005,0.01,0.02,0.05,0.1,0.2,0.5, 1, 2, 3]';
    sweepvals_Loss_SH_dBcm  =  [10,20,30,40,50,100,200,400,600,800,1000,2000,3000,30000,100000]';
    sweepvals_Loss_SH_dBcm  = [10,20,30,40,50,100,200,400,600,800,1000,2000,3000,5000,10000,20000,30000,100000]';
    FOM = zeros(length(sweepvals_delta_lam),length(sweepvals_Loss_SH_dBcm));
    for ksweep1 = 1:length(sweepvals_delta_lam)
        modn.delta_lam = sweepvals_delta_lam(ksweep1);
        for ksweep2 = 1:length(sweepvals_Loss_SH_dBcm)
            modn.Loss_SH_dBcm = sweepvals_Loss_SH_dBcm(ksweep2);
            [tmpmodn,PscattRatio,P3guided_max, P3guided_end, P3b_scatt_end, P3b_scatt_max] = integrate_z(c,modn,0);
            FOM(ksweep1,ksweep2) = PscattRatio;
            P3ss(ksweep1,ksweep2) = P3guided_max;
            %P3scatt(ksweep1,ksweep2) = P3b_scatt_end;
            P3scatt(ksweep1,ksweep2) = P3b_scatt_max;
        end
    end
    
    %% P3_guided (end)
   figure; 
    pcolor(sweepvals_Loss_SH_dBcm,sweepvals_delta_lam,log10(P3ss));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('detuning [nm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title(['P3(guided): max value (0 - ',num2str(modn.L_waveguide_um),')']);
   
    %% P3_scatt (end)
   figure; 
    pcolor(sweepvals_Loss_SH_dBcm,sweepvals_delta_lam,log10(P3scatt));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('detuning [nm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title('P3(scatt): total over full length')
    %caxis([-10,-3]);
    
   figure; 
    surf(sweepvals_Loss_SH_dBcm,sweepvals_delta_lam,log10(P3scatt));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('detuning [nm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title('P3(scatt): total over full length')
    %caxis([-10,-3]);

    
    
    %% FOM
    figure; 
    pcolor(sweepvals_Loss_SH_dBcm,sweepvals_delta_lam,log10(FOM));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('detuning [nm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title('Fraction of Generated Photons that are out-coupled')
    
   
    %linewidth calculation (threshold)
    for ksweep2 = 1:length(sweepvals_Loss_SH_dBcm)
        Y(:,ksweep2) = FOM(:,ksweep2)/ max(FOM(:,ksweep2));
        LW(ksweep2)= interp1(Y(:,ksweep2),sweepvals_delta_lam,0.5);
        legnames{ksweep2} = [num2str(sweepvals_Loss_SH_dBcm(ksweep2)),'dB/cm'];
    end
    
    figure; loglog(sweepvals_Loss_SH_dBcm,LW);
    xlabel('Signal (SH) out-coupling rate [dB/cm]');
    ylabel('detuning [nm] ~  HWHM ');
    
    figure;
    loglog(sweepvals_delta_lam,Y');  xlabel('detuning [nm]');  ylabel('normalized FOM')
    legend(legnames);
    
    figure;
    semilogx(sweepvals_delta_lam,FOM');  xlabel('detuning [nm]');  ylabel('Fraction of Generated Photons that are out-coupled');
    legend(legnames);
    
    figure;
    semilogx(sweepvals_delta_lam,Y');  xlabel('detuning [nm]');  ylabel('normalized FOM');
    legend(legnames);
       
    disp(sweepvals_Loss_SH_dBcm);
    disp(LW');
    
    close all;
    iselp = [1,4,5,8,11,15];
    figure; plot(sweepvals_delta_lam,P3scatt(:,iselp)');  xlabel('detuning [nm]');  ylabel('P3Scatt');
    legend(legnames{iselp}); grid on; xlim([0.005,1]); 
    figure; semilogx(sweepvals_delta_lam,P3scatt(:,iselp)');  xlabel('detuning [nm]');  ylabel('P3Scatt');
    legend(legnames{iselp}); grid on; xlim([0.005,1]);
    figure; semilogy(sweepvals_delta_lam,P3scatt(:,iselp)');  xlabel('detuning [nm]');  ylabel('P3Scatt');
    legend(legnames{iselp}); grid on; xlim([0.005,1]); ylim([1e-7,1e-4])   
    
    figure; plot(sweepvals_delta_lam,P3ss(:,iselp)');  xlabel('detuning [nm]'); ylabel('P3guided max(0...L)');
    legend(legnames{iselp}); grid on; xlim([0.005,1]); 
    figure; semilogx(sweepvals_delta_lam,P3ss(:,iselp)');  xlabel('detuning [nm]');  ylabel('P3guided max(0...L)');
    legend(legnames{iselp}); grid on; xlim([0.005,1]);  ylim([0,1.5e-4]) 
    figure; semilogy(sweepvals_delta_lam,P3ss(:,iselp)');  xlabel('detuning [nm]');  ylabel('P3guided max(0...L)');
    legend(legnames{iselp}); grid on; xlim([0.005,1]); ylim([1e-7,1e-4])   
    
    
    
    
    
    
    
    %%
    %-----------------------------------------------------------------------------------------------------------------
    %sweep: ABSORPTION3 & ALPHA3_scatter
    %-----------------------------------------------------------------------------------------------------------------
    clear FOM P3ss P3scatt;
    modn.delta_lam            = 0.5;      %this is delta from ideal PhaseMatch WL.  
    sweepvals_ABS_SH_dBcm     = [0.1,1,2,10,20,30,40,50,100,200,500,1000]';
    sweepvals_Loss_SH_dBcm  = [10,20,30,40,50,100,200,400,600,800,1000,2000,3000,30000,100000]';
    FOM = zeros(length(sweepvals_ABS_SH_dBcm),length(sweepvals_Loss_SH_dBcm));
    for ksweep1 = 1:length(sweepvals_ABS_SH_dBcm)
        modn.Abs_SH_dBcm = sweepvals_ABS_SH_dBcm(ksweep1);
        for ksweep2 = 1:length(sweepvals_Loss_SH_dBcm)
            modn.Loss_SH_dBcm = sweepvals_Loss_SH_dBcm(ksweep2);
            [tmpmodn,PscattRatio,P3guided_max, P3guided_end, P3b_scatt_end] = integrate_z(c,modn,0);
            FOM(ksweep1,ksweep2) = PscattRatio;
            P3guided(ksweep1,ksweep2) = P3guided_max;
            P3scatt(ksweep1,ksweep2) = P3b_scatt_end;
        end
    end
    
    
    %% P3_guided (end)
   figure; 
    pcolor(sweepvals_Loss_SH_dBcm,sweepvals_ABS_SH_dBcm,log10(P3guided));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('Absorp Loss [dB/cm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title(['P3(guided): max value (0 - ',num2str(modn.L_waveguide_um),')']);
   
   caxis([-9,-5]);
    
    %% P3_scatt (end)
   figure; 
    pcolor(sweepvals_Loss_SH_dBcm,sweepvals_ABS_SH_dBcm,log10(P3scatt));
    shading interp; colorbar;
    ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log');
    ylabel('Absorp Loss [dB/cm]'); xlabel('Signal (SH) out-coupling rate [dB/cm]');
    title('P3(scatt): total over full length')
   caxis([-9,-5]);
     
    %linewidth calculation (threshold)
    for ksweep2 = 1:length(sweepvals_Loss_SH_dBcm)
        Yscatt(:,ksweep2) = P3scatt(:,ksweep2)./ max(P3scatt(:,ksweep2));
        LWs(ksweep2)= interp1(Yscatt(:,ksweep2),sweepvals_ABS_SH_dBcm,0.9);
        
        Yg(:,ksweep2) = P3guided(:,ksweep2)./ max(P3guided(:,ksweep2));
        LWg(ksweep2)= interp1(Yg(:,ksweep2),sweepvals_ABS_SH_dBcm,0.9);
        legnames{ksweep2} = [num2str(sweepvals_Loss_SH_dBcm(ksweep2)),'dB/cm'];
    end
   
    [dummy,isel_g] = max(P3guided(1,:));
    [dummy,isel_scatt] = max(P3scatt(1,:));
    figure; semilogx(sweepvals_ABS_SH_dBcm,[Yg(:,isel_g),Yscatt(:,isel_scatt)]);  xlabel('ABSORPTION Loss [dB/cm]');  ylabel('Relative Power')
    xlabel('ABS [dB/cm]');   ylabel('Relative Power');  grid on;
    legend(['Guided: (at scatt=',num2str(sweepvals_Loss_SH_dBcm(isel_g)),'dB/cm)'],...
           ['Radiated: (at scatt=',num2str(sweepvals_Loss_SH_dBcm(isel_scatt)),'dB/cm)'] );  
     title([ 'Tolerance of ABS Loss (for optimal scatter rate).  detune=',num2str(modn.delta_lam),'nm']);
     
%     figure; loglog(sweepvals_Loss_SH_dBcm,[LWg',LWs']);
%     xlabel('Signal (SH) out-coupling rate [dB/cm]');
%     ylabel('linewidth [nm]');


end

    



% function dYdz = Yvec_primeB(z,Y, a, b, dk)
%     %     d2y/dz2 + a*dy/dz + b*exp(i*dk*z) = 0;    %ORIGINAL EQUATION
%     %     dA1/dz = -alpha0*A1;                      %ORIGINAL EQUATION
%     
%     %     dy1/dz = y2;                              %definition of y2
%     %     dy2/dz = - a*y2 - b*exp(i*dk*z);          %second order equation converted to first order    
%     dYdz = [Y(2);  -a*Y(2) - b*exp(i*dk*z) ];        
% end
function dYdz = Yvec_prime(z,Y, a, b, dk)
    %     d2y/dz2 + a*dy/dz + b*exp(i*dk*z) = 0;    %ORIGINAL EQUATION
    %     dy1/dz = y2;                              %definition of y2
    %     dy2/dz = - a*y2 - b*exp(i*dk*z);          %second order equation converted to first order    
   dYdz = [Y(2);  -a*Y(2) - b*exp(i*dk*z) ];        
end


 
function [modn,PscattRatio, P3guided_max, P3guided_end, P3scatt_end, P3scatt_max ] = integrate_z(c,modn,doplot)
  
    %parameters calc
    modn.lam2            = modn.lam1 / 2;                           %[nm] 
    modn.w_1             = 2*pi* (c.c0) / (  modn.lam1 * 1e-9);    %[Hz]
    modn.w_2             = 2*pi* (c.c0) / (  modn.lam2 * 1e-9);    %[Hz]  
    modn.prefactor       = 2*modn.w_2^2*(modn.d_xx*1e-12)^2 / (c.eps0*modn.n1^2*modn.n2*c.c0^3);    %[um]
    modn.A_waveguide_um2 = modn.x_waveguide_um * modn.y_waveguide_um;                               %[um^2]
    modn.g               = sqrt( modn.overlap_frac^2 / modn.A_waveguide_um2 * modn.prefactor );     %[um-1]
    
    % Detuning (Center of Pump Line) +  Offset from Linewidth (modn.linewidth_nm/2) * 2;
    delta_n     = modn.delta_lam * modn.dispersion_slope_diff;    %unitless
    delta_k_cm  = 4*pi/(modn.lam1*1e-7) * delta_n;           %[cm^-1]
    modn.alpha0_cm = modn.Loss_pump_dBcm / 4.3429;           % pump: alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
    modn.alpha3_cm = modn.Loss_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
    modn.abs3_cm   = modn.Abs_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34

    a0_um      = modn.alpha0_cm  * (1e-4); 
    a3_um      = modn.alpha3_cm  * (1e-4); 
    dk_um      = delta_k_cm     * (1e-4);               %[um-1]
    g          = modn.g;                                    %[um-1]
    abs3_um    = modn.abs3_cm  * (1e-4);
    
    %numeric integration over 0...L
    dz_um = c.dz_um;                                        %stepsize [um]
	zvec = [0:dz_um:modn.L_waveguide_um]';
	%initialize
    A1 = zeros(size(zvec));
    A3 = zeros(size(zvec));
    A3scatt = zeros(size(zvec));
    P3_scatt_z = zeros(size(zvec));
    P3b_scatt_z = zeros(size(zvec));
        
    %initial conditions 
    phase1 = 0;  phase3 = 0;
    A1(1) = sqrt(modn.Pp)*exp(i*phase1);                 %units:[W^0.5]
    A3(1) = 0;                                           %[W^0.5]
    A3scatt(1) = 0;                                      %[W^0.5]
    A3_nophase(1) = 0;                                   %[W^0.5]
    P3_nophase(1) = 0;                                   %[W^0.5]
    P3_scatt_z(1) = 0;
    P3b_scatt_z(1) = 0;
    dP3_loss = 0; 
    dP3_gain=0;

    %Euler Setup
    n3 = modn.modeindex_crossing;
    k3_um = 2*pi/(modn.lam2*1e-3)*n3;
    a = 2*i*k3_um;                                  %[um-1]
    euler_k.y  = 0;                                 %initial condition:  y(0) = 0;
    euler_k.y1 = 0;                                 %initial condition:  y'(0) = 0;
    Y(1:3,1) = 0;

%     close all;
%     zspan   = [0,1000];
%     Y0      = [0;0];
%     b       = i*g*abs(A1(1,1))^2;
%     a       = 2*i*k3_um;            %[um-1]
%     [zvec,Y] = ode45(@(z,Y) Yvec_prime(z,Y,a,b,dk_um),zspan,Y0);
%     [zvec,YY] = ode45(@(z,YY) Yvec_primeB(z,YY,a,g,dk_um),zspan,Y0);
%      
    
    close all;
    for kz = 1:length(zvec)-1            
        
        %PUMP loss
        dA1_loss(kz,1) = -0.5*a0_um*A1(kz,1) * dz_um;            %scattering (at SH)
        
        %SH generation
        dA1genn(kz,1) =  - g * (abs(A1(kz,1))^2) * dz_um;
        dA3genn(kz,1) = i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz)) * dz_um;
        
        %SH Loss
        dA3loss(kz,1) = -0.5*a3_um*A3(kz,1) * dz_um;            %scattering (at SH)
        dA3abs(kz,1) = -0.5*abs3_um*A3(kz,1) * dz_um;           %absorption (at SH)
        
        %Scatter Collection
        dP3b_scatt = (1-exp(-a3_um*dz_um))*(A3(kz,1)*conj(A3(kz,1)));       %units[W]
        %         dP3_scatt(kz,1) = dA3loss(kz,1)*conj(dA3loss(kz,1));            % ???? WRONG?? not sure magnitude is correct ...
    	
        %INTEGRATE
        A3_nophase(kz+1,1) = A3_nophase(kz,1) + abs(dA3genn(kz,1));              %add up all photons generated, **ignoring phase**
        A1(kz+1,1) = A1(kz) + dA1_loss(kz,1) + dA1genn(kz,1);                    %pump: guided     (units^0.5])
        A3(kz+1,1) = A3(kz,1) + dA3genn(kz,1) + dA3loss(kz,1) + dA3abs(kz,1)  ;  %SH: guided     (units^0.5])                        
        P3b_scatt_z(kz+1,1) = P3b_scatt_z(kz,1) + dP3b_scatt;                    %SH:  scattered   (units[W])
        %         P3_nophase(kz+1,1) = P3_nophase(kz,1) + abs(dA3genn(kz,1))^2;     %add up all photons generated, **ignoring phase**
        
%         %%2nd Order Euler Method %%%%%%%%%%%%%%
%         b = i*g*abs(A1(kz,1))^2 ;
%          dYdz(kz,1) = 
%          Y(kz+1,1) = Y(kz,1) + h*dYdz(kz,1);
%          Y(kz+1,2) = Y(kz,2) + h*dYdz(kz,2);         
        
%         %%2nd Order Euler Method %%%%%%%%%%%%%%
%         b = i*g*abs(A1(kz,1))^2 ;
%         euler.A3(kz,1) = euler_k.y;
%         euler_k.y2 =  +b*exp(-i*dk_um*zvec(kz)) - a* euler_k.y1; 
%         euler_kplus1.y  = euler_k.y  + euler_k.y1 * dz_um;
%         euler_kplus1.y1 = euler_k.y1 + euler_k.y2 * dz_um;
%         euler_k = euler_kplus1;                                     %step forward
% 
%         
%         d2y_dz2 =  - a*Y(2,kz)  + i*g*abs(A1(kz,1))^2*exp(-i*dk_um*zvec(kz)) + dA3loss(kz,1) + dA3abs(kz,1);
%         Y(3,kz) = d2y_dz2;
%         Y(1,kz+1) = Y(1,kz) + Y(2,kz) * dz_um;
%         Y(2,kz+1) = Y(2,kz) + Y(3,kz) * dz_um;
        
     end %zloop
%     euler.P3 = euler_k.y .* conj(euler_k.y);
    P1 = A1.*conj(A1); 
    P3 = A3.*conj(A3);
    P3totalgenerated = (A3_nophase).^2;
    PscattRatio = P3b_scatt_z(end)/P3totalgenerated(end);
    P3guided_end = P3(end);
    P3guided_max = max(P3);
    P3scatt_end = P3b_scatt_z(end);
    P3scatt_max = max(P3b_scatt_z);
    
    if doplot==1
        figure; 
%         subplot(4,1,1); plotyy(zvec,[P3],zvec,[P3totalgenerated,P3b_scatt_z]); ylabel('Power [W]');  legend('P3_g_u_i_d_e_d','P3_t_o_t_a_l_g_e_n_e_r_a_t_e_d','P3b_s_c_a_t_t')
        subplot(4,1,1); plotyy(zvec,[P1,P3,P3totalgenerated],zvec,[P3b_scatt_z]); ylabel('Power [W]');  legend('P3_g_u_i_d_e_d','P3_t_o_t_a_l_g_e_n_e_r_a_t_e_d','P3b_s_c_a_t_t')
        title({['delta k =',num2str(round(delta_k_cm,1)),'[cm-1]  alpha3=',num2str(round(modn.alpha3_cm)),'[cm-1]'];...
            ['detuning: ',num2str(modn.delta_lam),'nm  P3scatt(L)/P3totgen(L)=',num2str(P3b_scatt_z(end)/P3totalgenerated(end))]});    
        subplot(4,1,2);  semilogy(zvec,[P3./P3totalgenerated,P3b_scatt_z./P3totalgenerated]); ylabel('Fraction of Totalgenerated')
        legend('P3_g_u_i_d_e_d ~ SINC(z)','P3b_s_c_a_t_t')
        subplot(4,1,3); plotyy(zvec,[abs(A3)],zvec(1:end-1),[abs(dA3genn),abs(dA3loss)]); ylabel('amplitiude');
        legend('A3','dA3gen','dA3loss');
        subplot(4,1,4); plot(zvec(1:end-1),[rad2deg(unwrap(angle(A3(1:end-1)))), rad2deg(unwrap(angle(dA3genn))),rad2deg(unwrap(angle(dA3loss)))]); 
        ylabel('phase[deg]');
        xlabel('z[um]');
        
      
        figure; 
        semilogy(zvec,[P3./P3totalgenerated,P3b_scatt_z./P3totalgenerated]); ylabel('Fraction of Totalgenerated')
        legend('P3_g_u_i_d_e_d ~ SINC(z)','P3b_s_c_a_t_t');  ylabel('P3b scatt_z / P3totalgenerated');        
    end
    
    
    
end
    


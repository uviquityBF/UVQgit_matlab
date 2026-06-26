function calc_numeric_zerolinewidth_TEST()

    %define inputs
        c.c0                        = 3e8; %[m/s]
        c.eps0                      = 8.85E-12;  %[C/(V.m)]
        c.dz_um                     = 0.25;         %[um]  stepsize for numeric integration 
        
        modn.dispersion_slope_diff  = 0.009;
        modn.name                   = '20 dB/cm';
        modn.Pp                     = 0.100000000000000;
        modn.linewidth_nm           = 1.000000000000000e-05;
        modn.d_xx                   = 4.700000000000000;
        modn.lam1                   = 450;
        modn.n1                     = 2.150000000000000;
        modn.n2                     = 2.600000000000000;
        modn.overlap_frac           = 0.030000000000000;
        modn.x_waveguide_um         = 0.336500000000000;
        modn.y_waveguide_um         = 0.381000000000000;
        modn.dispersion_slope_diff  = 0.008950439850000;
        
        modn.L_waveguide_um         = 10000;
        modn.Loss_pump_dBcm         =  20;
        modn.Loss_SH_dBcm           = 200;
        
        modn.Pump_Detune_nm         = 0;
        modn.lam2                   = 225;
        modn.w_1                    = 1.948274513854135e+15;
        modn.w_2                    = 3.222146311374147e+15;
        modn.prefactor              = 1.597198607427036e-07;
        modn.A_waveguide_um2        = 0.128206500000000;
        modn.g                      = 3.348464462025136e-05;

    %call and test function
        delta_lam_linewidth_nm = 0;  %Offset from Center of Linewidth (modn.linewidth_nm/2) * 2
        modn = calc_numeric_zerolinewidth(c,modn,delta_lam_linewidth_nm);

end


function  modn = calc_numeric_zerolinewidth(c,modn,delta_lam_linewidth_nm)

       %run NUMERIC Integation  ..................................
        
        %parameters define:
        delta_lam = modn.Pump_Detune_nm + delta_lam_linewidth_nm;     %this is delta from ideal PhaseMatch WL.  it is sum of:
                                                                 % Detuning (Center of Pump Line) +  Offset from Center of Linewidth (modn.linewidth_nm/2) * 2;
        delta_n     = delta_lam * modn.dispersion_slope_diff;    %unitless
        delta_k_cm  = 4*pi/(modn.lam1*1e-7) * delta_n;           %[cm^-1]
        modn.alpha0_cm = modn.Loss_pump_dBcm / 4.3429;           % pump: alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
        modn.alpha3_cm = modn.Loss_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
        dk_um      = delta_k_cm     * (1e-4);
        a0_um      = modn.alpha0_cm  * (1e-4); 
        a3_um      = modn.alpha3_cm  * (1e-4); 
        g          = modn.g;

        %numeric integration over 0...L
        dz_um = c.dz_um;                                        %stepsize [um]
        zvec = [0:dz_um:modn.L_waveguide_um]';              
        %initialize
        A1 = zeros(size(zvec));
        A3 = zeros(size(zvec));
        A3scatt = zeros(size(zvec));
        P3b_scatt_z = zeros(size(zvec));
        
        %initial conditions 
        phase1 = 0;  phase3 = 0;
        A1(1) = sqrt(modn.Pp)*exp(i*phase1);                 %units:[W^0.5]
        A3(1) = 0;                                           %[W^0.5]
        A3scatt(1) = 0;                                      %[W^0.5]
        P3b_scatt_z(1) = 0;
        P1b_scatt_z(1) = 0;
        modtmp=modn;
        for kz = 1:length(zvec)-1      

            %diffeq_1:    dA3/dz = -0.5a3*A3  + i*g*A1^2*exp(-i*dk*dz)
            dA3 = ( -0.5*a3_um*A3(kz,1)   +    i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz)) ) * dz_um;
            A3(kz+1,1) = A3(kz,1) + dA3;

            %diffeq_2:    dA1/dz = -0.5a0*A1 ;
            dA1 = -0.5*a0_um*A1(kz,1) * dz_um;
            A1(kz+1,1) = A1(kz,1) + dA1;                                        %units:[W^0.5]

            %diffeq_3:    dA3scatt/dz = +0.5a3*A3*dz;
            dA3scatt = +0.5*a3_um*A3(kz,1) * dz_um;
            A3scatt(kz+1,1) = A3scatt(kz,1) + dA3scatt;                         %units:[W^0.5]
               
            %diffeq_3B:    dP3scatt/dz = (1-exp(a3dz))*P3;
            dP3b_scatt = (1-exp(-a3_um*dz_um))*(A3(kz,1)*conj(A3(kz,1)));       %units[W]
            P3b_scatt_z(kz+1,1) = P3b_scatt_z(kz,1) + dP3b_scatt;               %units[W] 

            %diffeq_1B:    dP1scatt/dz = (1-exp(a1dz))*P1;
            dP1b_scatt = (1-exp(-a0_um*dz_um))*(A1(kz,1)*conj(A1(kz,1)));       %units[W]
            P1b_scatt_z(kz+1,1) = P1b_scatt_z(kz,1) + dP1b_scatt;               %units[W] 
            
        end %Integrate_over_Z
        
        % Amplitude --> Power
        P1z                 = A1.*conj(A1);                                   %[W]
        P3z                 = A3.*conj(A3);                                   %[W]
        P3scatt_z           = A3scatt.*conj(A3scatt);                         %[W]

        %for energy conservation check
        Pp_z                = modn.Pp * ones(size(zvec));   
        Psum_z              = P1z + P3z + P1b_scatt_z + P3scatt_z;            %all outputs (using P3scatt)
        Psum_b_z             = P1z + P3z + P1b_scatt_z + P3b_scatt_z;        %all outputs (using P3b_scatt)
        
        % PLOT: RESULTS
            figure; plotyy(zvec, P1z, zvec, [P3z, P3scatt_z, P3b_scatt_z] ); ylabel('power [W]'); xlabel('z [um]');
            legend('P1','P3','P3scatt','P3bscatt(Wei)');    title({'Numeric Results (only)',['delta WL = ',num2str(delta_lam),'[nm]']});
            figure; semilogy(zvec, [P1z, P3z, P3scatt_z, P3b_scatt_z] ); ylabel('power [W]'); xlabel('z [um]');
            legend('P1','P3','P3scatt','P3bscatt(Wei)');    title({'Numeric Results (only)',['delta WL = ',num2str(delta_lam),'[nm]']});
             
        %energy conservation check -- should never exceed Pp
        figure; semilogy(zvec,[Pp_z,Psum_z,Psum_b_z]); ylabel('power [W]'); xlabel('z [um]');
        legend('P_p_u_m_p','P_S_U_M(all outputs)','P_S_U_M b (all outputs)'); 
        title({'Energy Conservation Check';'Expect slight error because Pump depletion NOT modeled'});
        ylim([0.09999,0.1001])
        
end  %NUMERIC CALC
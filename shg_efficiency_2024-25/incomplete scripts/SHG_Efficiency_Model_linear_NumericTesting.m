%  SHG Efficiency Model(s)
%  based : Boyd 2.2.19 (with Area = Area_waveguide)
%  7/15/2025 -- NUMERIC INTEGRATION
%  TEMP -- COMPARISON OF NUMERIC v. ANALYTIC

close all;
clear;
c.c0                = 3e8; %[m/s]
c.eps0              = 8.85E-12;  %[C/(V.m)]

% dirs.h = 'C:\Users\brent\MATLAB\matDat';
dirs.h = 'G:\My Drive\Analyses(BF)\SHG Ideation\SHG_Efficiency_Models&Formulas';
inputfile =  'SHGEfficModel(matlab)_INPUT.xlsx';
outputfile = 'temp.csv';

%% Input Data  (separate "case" for each column)
INs = get_INPUT(dirs,inputfile);


%% set up 
IN = INs{4};
mod = IN;       %copy input values to this model ('mod')
  


%% DEFINE SWEEPsover Selected Parameters
 
sweep1.fieldname  =  'L_waveguide_um';
sweep1.values      = [1000:1000:2000]';  %[W]
%%
 %LOOP OVER CASES of INPUT PARAMETERS (columns of input)
MR = []; 
do_plots=0;
for kc=3:length(INs)

    legnames{kc} = INs{kc}.name;
     
    %SWEEP Z  = LENGTH  ...........................
    for k = 1:length(sweep1.values)     
      %reset "IN" + change Parameter
        IN = INs{kc};  
        IN.(sweep1.fieldname)           =  sweep1.values(k,1)';
        mod = IN;

      %run ANALYTIC Model at ZERO LINEWIDTH
        sweep1.mod(k)                   =  calcAnalytic_ZeroLinewidth(c,IN,mod,do_plots);  
        
      %Record key results for ANALYTIC
        sweep1.a.Effic_ZeroLinewidth_noDetune_NoLoss(k,kc) ...
                                        = sweep1.mod(k).Efficiency_ZeroLinewidth_noDetune_NoLoss;
        sweep1.a.DeRate_guided(k,kc)    = sweep1.mod(k).DeRate_Wei;
        sweep1.a.Effic(k,kc)          	= sweep1.mod(k).Efficiency_ZeroLinewidth_noDetune_NoLoss  *  sweep1.mod(k).DeRate_Wei;
        sweep1.a.Psh(k,kc)              = sweep1.a.Effic(k,kc) * sweep1.mod(k).Pp;    
        
      %run NUMERIC Integation  ..................................
        modn = calc_prefactor(c,mod);      %modn  holds data for *NUMERIC* integration result   
        
        %parameters define:
        delta_lam   = modn.Pump_Detune_nm;                       %(modn.linewidth_nm/2) * 2;
        delta_n     = delta_lam * modn.dispersion_slope_diff;    %unitless
        delta_k_cm  = 4*pi/(modn.lam1*1e-7) * delta_n;           %[cm^-1]
        modn.alpha0_cm = modn.Loss_pump_dBcm / 4.3429;           % pump: alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
        modn.alpha3_cm = modn.Loss_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
        dk_um      = delta_k_cm     * (1e-4);
        a0_um      = modn.alpha0_cm  * (1e-4); 
        a3_um      = modn.alpha3_cm  * (1e-4); 
        g          = modn.g;

        %initial conditions 
        phase1 = 0;  phase3 = 0;
        A1(1) = sqrt(modn.Pp)*exp(i*phase1);                 %units:[W^0.5]
        A3(1) = 0;                                           %[W^0.5]
        A3scatt(1) = 0;                                      %[W^0.5]

        %numeric integration over 0...L
        dz_um = 0.1;                                        %stepsize [um]
        zvec = [0:dz_um:modn.L_waveguide_um]';
        modtmp=mod;
        for kz = 1:length(zvec)-1      

            %diffeq_1:    dA3/dz = -0.5a3*A3  + i*g*A1^2*exp(-i*dk*dz)
            dA3 = ( -0.5*a3_um*A3(kz,1)   +    i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*dz_um) ) * dz_um;
            A3(kz+1,1) = A3(kz,1) + dA3;

            %diffeq_2:    dA1/dz = -0.5a0*A1 ;
            dA1 = -0.5*a0_um*A1(kz,1) * dz_um;
            A1(kz+1,1) = A1(kz,1) + dA1;                    %units:[W^0.5]

            %diffeq_3:    dA3scatt/dz = +0.5a3*A3;
            dA3scatt = +0.5*a3_um*A3(kz,1) * dz_um;
            A3scatt(kz+1,1) = A3scatt(kz,1) + dA3scatt;          %units:[W^0.5]
               
            %analytic result:  def_Integral(0...L)
            L                       = zvec(kz);
            modtmp.L_waveguide_um   = L;
            modtmp = calc_prefactor(c,modtmp); 
            modtmp = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,modtmp);        %
            modtmp = PhaseMismatch_and_Loss_Derate(c, modtmp, modtmp.Pump_Detune_nm);   % analytic solution
            
            Eff_zlwnl(kz+1,1)               = modtmp.Efficiency_ZeroLinewidth_noDetune_NoLoss;
            Eff_guided(kz+1,1)              = modtmp.DeRate_Wei * Eff_zlwnl(kz+1,1);              
            Eff_scatt(kz+1,1)               = modtmp.DeRate_Scatt * Eff_zlwnl(kz+1,1);               
            P3z_analytic(kz+1,1)            = modtmp.Pp * Eff_guided(kz+1,1);           %[W]
            PscattONLY_analytic(kz+1,1)     = modtmp.Pp * Eff_scatt(kz+1,1);            %[W]
            P3plusScatt_analytic(kz+1,1)    = modtmp.Pp * Eff_scatt(kz+1,1);            %[W]
 
        end %Integrate_over_Z
        
        % Amplitude --> Power
        P1z                 = abs(A1).^2;                                   %[W]
        P3z                 = abs(A3).^2;                                   %[W]
        P3scatt_z           = abs(A3scatt).^2;                              %[W]
        
        % COMPARE:  NUMERICAL v. ANALYTICAL
        figure; plotyy(zvec, P1z, zvec, [P3z, P3scatt_z] ); ylabel('power [W]'); xlabel('z [um]');
        legend('P1','P3','P3scatt');    
        
        size(zvec)
        size([P3z, P3z_analytic, P3scatt_z, PscattONLY_analytic])
        figure;semilogy(zvec, [P3z, P3z_analytic, P3scatt_z, PscattONLY_analytic]);  ylabel('power [W]'); xlabel('z [um]');
        legend('guided: numeric','guided: analytic','scatt(Only): numeric','scatt(Only): analytic');  
        

        %Value at z=L
        modn.P1             = P1z(end);
        modn.Psh_guided     = P3z(end);
        modn.Psh_wscatt 	= P3z(end)+P3scatt_z(end);
        
        %Record key results for NUMERICAL
        sweep1.n.Psh_guided(k,kc)       = modn.Psh_guided;    
        sweep1.n.Psh_wscatt(k,kc)       = modn.Psh_wscatt;    
        sweep1.n.Effic_guided(k,kc)     = (modn.Psh_guided) / modn.Pp;    
        sweep1.n.Effic_wscatt(k,kc)    	= (modn.Psh_guided + modn.Psh_wscatt) / modn.Pp;    
        

           
        
      % !!compare!!
        disp([sweep1.a.Psh(k,kc), sweep1.n.Psh_guided(k,kc)]);
      % IMPORTANT!!  CLEAR INEGRALS, else Vector(end) will not change on each sweep
        clear A1 A3 A3scatt P3z P1z P3scatt_z ;  
    end %end of Sweep    
     
end  %loop over cases


%PLOT:  shg power  v Length
hf=figure; semilogy(sweep1.values, sweep1.a.Psh*1e6,'Linewidth',2); grid on;
hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
semilogy(sweep1.values,  sweep1.n.Psh_guided*1e6,'*','MarkerSize',4);  
hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
semilogy(sweep1.values,  sweep1.n.Psh_wscatt*1e6,'--','Linewidth',4);  
xlabel(replace(sweep1.fieldname,'_',' ')); ylabel('P_S_H (avg) [uW] ');
legend([legnames]);  title(replace(sweep1.fieldname,'_',' '));


%Write Results for All Cases
csvwrite([dirs.h,'/',outputfile], MR);
disp(['Wrote output to:', outputfile ]);



%%  Calc All
function mod = calcAnalytic_ZeroLinewidth(c,IN,mod,do_plots)
    mod = IN;
    mod = calc_prefactor(c,mod); 
    mod = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,mod);      %
    mod = PhaseMismatch_and_Loss_Derate(c, mod, mod.Pump_Detune_nm);    %analytic solution
    % -------------------
%     
%     mod.Lcoh_nm = -1;
%     mod.DeRate_Overlap_guided = 1;
%     mod.DeRate_Overlap_ScattOnly = -1;
%     mod.Efficiency_ScattOnly = -1;
%     mod.DerateUsed = 'DeRate_Wei';
%     mod.Efficiency = mod.GAMMA_pctPerWatt * mod.Pp^2 * mod.DeRate_Wei;    
%     mod.Efficiency_ScatterAndGuided = -1;
%     mod.Psh         = mod.Efficiency * mod.Pp;
%     % -------------------
end


%%  Prefactor Calc
function mod = calc_prefactor(c,mod) 
    mod.lam2            = mod.lam1 / 2;
    mod.w_1             = 2*pi* (c.c0) / (  mod.lam1 * 1e-9);    %[Hz]
    mod.w_2             = 2*pi* (c.c0) / (  mod.lam2 * 1e-9);    %[Hz]
    mod.prefactor       = 2*mod.w_2^2*(mod.d_xx*1e-12)^2 / (c.eps0*mod.n1^2*mod.n2*c.c0^3);
    mod.A_waveguide_um2  = mod.x_waveguide_um * mod.y_waveguide_um;     %[um^2]
    mod.g               = sqrt( mod.overlap_frac^2 / mod.A_waveguide_um2 * mod.prefactor );        % UPDATE THIS **
end
%%
function mod = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,mod)
    mod.A_waveguide_um2  = mod.x_waveguide_um * mod.y_waveguide_um;  %[um^2]
    mod.GAMMA_pctPerWatt = mod.prefactor * mod.overlap_frac^2 ...
                            * (mod.L_waveguide_um^2 / mod.A_waveguide_um2);  
    mod.Efficiency_ZeroLinewidth_noDetune_NoLoss    = mod.GAMMA_pctPerWatt * mod.Pp;
end

%%
function mod =  PhaseMismatch_and_Loss_Derate(c, mod, delta_lam_nm)
    delta_lam   = delta_lam_nm;                             %(mod.linewidth_nm/2) * 2;
    delta_n     = delta_lam * mod.dispersion_slope_diff;    %unitless
    delta_k_cm  = 4*pi/(mod.lam1*1e-7) * delta_n;           %[cm^-1]
    mod.alpha0_cm = mod.Loss_pump_dBcm / 4.3429;           % pump: alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
    mod.alpha3_cm = mod.Loss_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
    
    dkL  = delta_k_cm   * (mod.L_waveguide_um * 1e-4);
    a3L   = mod.alpha3_cm * (mod.L_waveguide_um * 1e-4); 
    a0L   = mod.alpha0_cm * (mod.L_waveguide_um * 1e-4); 
    daL   = a3L/2 - a0L;    
    
    %*** GUIDED LIGHT (Only) ****
    if dkL==0   
        mod.DeRateAppx0     = (exp(-a3L));
        mod.DeRateAppx      = (exp(-a3L));
        mod.DeRateAppx2     = (exp(-a3L));
        mod.DeRate_BF       = (exp(-a3L)) *2* ( cosh(a3L) - 1 ) / ( (a3L)^2 ) ;
        if daL==0
            mod.DeRate_Wei  = (exp(-a3L))*sinc(dkL/2/pi)^2;
        else
            mod.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
                                          / ( (dkL)^2 + (daL)^2 ) ;
        end
        %following 4 variables are here just to put them into struct
        mod.DeRate_Scatt            = []; 
        mod.DeRate_Scatt_indefInteg = [];
        mod.DeRate_Wei_IdealSHLoss  = [];
        mod.DeRate_Scatt_alt1       = [];
        if a0L==0 
            mod.DeRate_noPMreq  = 1;       
        else
            %original formula (brent's "incorrect" version):
            mod.DeRate_noPMreq  = (exp(-a0L)-1)^2 / ( (a0L)^2 ) ;       
        end      
        
    else %dkL>0
        % crude estimate based on SINC(x) only ....
        mod.DeRateAppx0 = (exp(-a3L)) *1* ( sinc(dkL/2) )^2 ;
        mod.DeRateAppx  = (exp(-a3L)) *1* ( sin(dkL/2) )^2 / (dkL/2)^2  ;
        mod.DeRateAppx2 = (exp(-a3L)) *1* ( sin(dkL/2) /  (dkL/2)  )^2  ;   
        
        %original formula (brent's "incorrect" version):
        mod.DeRate_BF     = (exp(-a3L)) *2* ( cosh(a3L) - cos(dkL) ) / ( (dkL)^2 + (a3L)^2 ) ;

        %formula from Wei:
        mod.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
                                          / ( (dkL)^2 + (daL)^2 ) ;   
        
        mod.DeRate_noPMreq = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*1 ) ...
                                          / ( (daL)^2 ) ;                                          
    end
      
    %*** SCATTERED LIGHT (Only) ***
    mod.DeRate_Scatt    = ((0.5*a3L)^2) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
                            % this results from using *DEFINITE*
                            % integration in the derivation
                            
    mod.DeRate_Scatt_indefInteg    ...
                            = ((0.5*a3L)^2) * ( exp(-2*a0L) ) ...
                            / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
                            % note: no sinc() function here -- (phasematch relaxed?)
                            % this results from using *INdefinite*
                            % integration in the derivation                  
                        
    mod.DeRate_Wei_IdealSHLoss  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / ( (dkL)^2 + (a0L)^2 );
                            %this sets alpha_SH Loss = 0 
                            % .... but still has phase matching requirement
            
    mod.DeRate_Scatt_alt1  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / ( (dkL)^2 + (a0L)^2 )   -  mod.DeRate_Wei;
                            % "Wei_Ideal  - Wei"
                        
  
end



%% get INPUT 
function INs = get_INPUT(dirs,inputfile)
    sheetname = 'Linear';
    %get input
    [tmp.num,tmp.txt,tmp.raw]       =  xlsread([dirs.h,'/',inputfile],sheetname, 'd5:z20' );
    [tmp2.num,tmp2.txt,tmp2.raw]    =  xlsread([dirs.h,'/',inputfile],sheetname, 'd4:z4' );

    if size(tmp.num,2)>=1 && size(tmp.num,1)>=10

       for kc=1:size(tmp.num,2)
            INs{kc}.name             = replace(tmp2.txt{kc},'_',' ');
            INs{kc}.Pp_avg_mW        = tmp.num(1,kc);   %[W] external pump (into waveguide), average power
            INs{kc}.DF               = tmp.num(2,kc);   %duty factor of pulses in the pump (if any) -- set to 1.0 for CW
            INs{kc}.Pp               = INs{kc}.Pp_avg_mW / INs{kc}.DF * 1e-3;    %instantaneous power of pulses  (if any)
            INs{kc}.linewidth_nm     = tmp.num(3,kc);   % [nm]
            INs{kc}.d_xx             = tmp.num(4,kc);   %pm/V
            INs{kc}.lam1             = tmp.num(5,kc);   
            INs{kc}.n1               = tmp.num(6,kc);   
            INs{kc}.n2               = tmp.num(7,kc);   
            INs{kc}.overlap_frac     = tmp.num(8,kc);  %fraction (of 1)
            INs{kc}.x_waveguide_um   = tmp.num(9,kc);  
            INs{kc}.y_waveguide_um   = tmp.num(10,kc); 
            INs{kc}.L_waveguide_um   = tmp.num(11,kc); 
            INs{kc}.dispersion_slope_diff   = tmp.num(12,kc); 
            INs{kc}.Loss_pump_dBcm   = tmp.num(13,kc); 
            INs{kc}.Loss_SH_dBcm     = tmp.num(14,kc); 
            INs{kc}.Pump_Detune_nm   = tmp.num(15,kc); 
        end
        Ncases = length(INs);
        
    else
        error('Excel Input was not obtained... Returning Default INPUTs ')
    end
    
end


%% Generate Output Matrix
function M_out = ResultsMatrix_AddColumn(MR, mod)
   newres      = zeros(5,1);
   newres(1,1) = mod.Lcoh_nm*1e-3;
   newres(2,1) = mod.prefactor;
   newres(3,1) = mod.GAMMA_pctPerWatt; 
   newres(4,1) = mod.Efficiency_ZeroLinewidth_noDetune_NoLoss;
   newres(5,1) = mod.DeRate_Overlap_guided;
   newres(6,1) = mod.Efficiency;
   newres(7,1) = mod.Efficiency_ScatterAndGuided;

    
   [nr,nc] = size(MR);
   if nr==0
       M_out = newres;
   else
       if nr == length(newres)
           M_out = [MR,newres];
       else
           error('Results List does not match size of Matrix, MR');
       end
   end   
end






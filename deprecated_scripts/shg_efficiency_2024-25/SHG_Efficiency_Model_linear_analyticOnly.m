%  SHG Efficiency Model(s)
%  based : Boyd 2.2.19 (with Area = Area_waveguide)

function main()

close all;
clear;
c.c0                            = 3e8; %[m/s]
c.eps0                          = 8.85E-12;  %[C/(V.m)]
c.dz_um                         = 0.05;     %[um]  stepsize for numeric integration 
c.N_steps_overlap_integration   = 41;       %make this ODD, so that a point appears at 0

c.test_derate_functions         = 0;        %generate plots of Analytic Derate = f(detune, alpha0) 
c.do_compare_to_analytic        = 0;        %show every numeric integration (PowerSH v Z):  analytic v numeric
c.do_plots                      = 0;        %Overlap_Plots



% dirs.h = 'C:\Users\brent\MATLAB\matDat';
dirs.h                = 'G:\My Drive\Analyses(BF)\SHG Ideation\SHG_Efficiency_Models&Formulas';
inputfile             = 'SHGEfficModel(matlab)_INPUT.xlsx';
sheetname             = 'Linear';
outputfile            = 'temp.csv';


%% Input Data  (separate "case" for each column)
INs     = get_INPUT(dirs,inputfile,sheetname);
IN      = INs{1};
moda    = IN;       %copy input values to this model ('mod')

if mod(c.N_steps_overlap_integration,2)==0  %make this ODD, so that a point appears at 0
    c.N_steps_overlap_integration=c.N_steps_overlap_integration+1;
end

%%  Test DeRate Functions
if c.test_derate_functions==1   test_DeRate_Functions(c,moda);  end

%% SINGLE RUN:  LOOP OVER CASES of INPUT PARAMETERS (columns of input)

MR = []; 
for kc=1:length(INs)
    %update Case
    IN = INs{kc}; 
    legnames{kc} = INs{kc}.name;
    
    %Run this Case
    moda = calcAnalytic(c,IN);
    modn = calcNumeric(c,IN);
    
    %compare numeric & Analytic  (guided)
    disp({'Analytic,    Numeric'; num2str([moda.Psh_guided,  modn.Psh_guided])})
    
    %Extract Results of Interest
    MR = ResultsMatrix_AddColumn(MR, moda, modn); 
end

%Write Results for All Cases
csvwrite([dirs.h,'/',outputfile], MR);
disp(['Wrote output to:', outputfile ]);


%% SWEEPS over Selected Parameters
Nsweeps = 1;        %!! How many (of below) sweeps to do!!

sweeps{1}.fieldname  =  'L_waveguide_um';
sweeps{1}.values      = linspace(100,3000,20)'; %[100:100:10000]';  %[um]

sweeps{2}.fieldname  =  'Pump_Detune_nm';
sweeps{2}.values     =  linspace(0, 1,  20)';  %[nm]

sweeps{3}.fieldname  =  'None';%'linewidth_nm';
sweeps{3}.values      = [1e-3,0.05,0.1: 0.1 :1.5]';  %[nm]

sweeps{3}.fieldname  =  'Loss_SH_dBcm';
sweeps{3}.values      = linspace(50,1000, 20)'; %[0:10:1000]';  %[dB/cm]

sweeps{5}.fieldname  =  'None';%'dispersion_slope_diff';
sweeps{5}.values      = linspace(0.1*moda.dispersion_slope_diff,2*moda.dispersion_slope_diff,10)';  %[1/nm]

 
%% PARAMETER SWEEPS

for ks = 1:Nsweeps
ks

%Parameter Sweep 1 ...........................
if ~strcmp(sweeps{ks}.fieldname,'None')
    for kc=1:length(INs)  
        legnames{kc} = INs{kc}.name;
        for k = 1:length(sweeps{ks}.values)
          %reset "IN" + change Parameter
            IN = INs{kc};  
            IN.(sweeps{ks}.fieldname)            =  sweeps{ks}.values(k,1)';
          %run ANALYTIC MODEL
            sweeps{ks}.a.mod(k)                   =  calcAnalytic(c,IN);       
          %Record key results
            sweeps{ks}.a.Effic_ZeroLinewidth(k,kc)= sweeps{ks}.a.mod(k).Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss;
            sweeps{ks}.a.Lcoh_um(k,kc)            = sweeps{ks}.a.mod(k).Lcoh_um;
            sweeps{ks}.a.DeRate_guided(k,kc)      = sweeps{ks}.a.mod(k).DeRate_Overlap_guided;  
            sweeps{ks}.a.DeRate_scatt(k,kc)       = sweeps{ks}.a.mod(k).DeRate_Overlap_ScattONLY;
            sweeps{ks}.a.Effic_guided(k,kc)       = sweeps{ks}.a.mod(k).Efficiency_guided;    
            sweeps{ks}.a.Effic_ScattONLY(k,kc)    = sweeps{ks}.a.mod(k).Efficiency_ScattONLY;
            sweeps{ks}.a.Effic_ScatterAndGuided(k,kc) = sweeps{ks}.a.mod(k).Efficiency_ScatterAndGuided;
            sweeps{ks}.a.Psh_guided(k,kc)         = sweeps{ks}.a.mod(k).Efficiency_guided * sweeps{ks}.a.mod(k).Pp;    
            sweeps{ks}.a.Psh_ScattONLY(k,kc)      = sweeps{ks}.a.mod(k).Efficiency_ScattONLY * sweeps{ks}.a.mod(k).Pp;    
            sweeps{ks}.a.Psh_ScattAndGuided(k,kc) = sweeps{ks}.a.mod(k).Efficiency_ScatterAndGuided * sweeps{ks}.a.mod(k).Pp;         
          %run NUMERIC MODEL
            c.do_plots=0;
            sweeps{ks}.n.modn(k) = calcNumeric(c,IN);           
           %Record key results      
            sweeps{ks}.n.Psh_guided(k,kc)         = sweeps{ks}.n.modn(k).Psh_guided;    
            sweeps{ks}.n.Psh_ScattONLY(k,kc)      = sweeps{ks}.n.modn(k).Psh_ScattONLY;    
            sweeps{ks}.n.Psh_ScattAndGuided(k,kc) = sweeps{ks}.n.modn(k).Psh_ScatterAndGuided;          
            
          %compare numeric & Analytic
            disp({'Analytic,    Numeric'; num2str([ sweeps{ks}.a.Psh_guided,   sweeps{ks}.n.Psh_guided])})

        end %sweep paramater values
    end %loop over cases
    
    %PLOTS of SWEEP RESULTS
    makeplot_SHGPower(legnames, sweeps{ks}, 'a');
%     makeplot_DeRates(legnames, sweeps{ks}, 'a') ;
%     makeplot_Efficiencies(legnames,sweeps{ks},'a');
    makeplot_Scatt2Guided(legnames,sweeps{ks},'a');
    
    makeplot_SHGPower(legnames, sweeps{ks}, 'n');
    makeplot_Scatt2Guided(legnames,sweeps{ks},'n');
    
    makeplot_ScattONLY_Power_CompareAnalyticNumeric(legnames, sweeps{ks}) ;
    
end %if~="None"
  
end %SWEEPS loop over "ks"



end %-- END: MAIN ----------------------------------------------------------







%% Plotting Functions
%--------------- Plot: De-Rates(detune,loss)
function makeplot_DeRates(legnames, sweepX, an) 
    if strcmp(an,'n') type='numeric'; else type='analytic'; end
    hf=figure; semilogy(sweepX.values, sweepX.(an).DeRate_guided ,'Linewidth',2); grid on;
    hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
    semilogy(sweepX.values, sweepX.(an).DeRate_scatt,'--','Linewidth',2);  
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('DeRate (process window) Value');
    legend([legnames]);  title({['DERATE Values (',type,')'];replace(sweepX.fieldname,'_',' ')});
   % xlim([0,2000]);  
end
%--------------- PowerSH  
function makeplot_SHGPower(legnames, sweepX, an) 
    if strcmp(an,'n') type='numeric'; else type='analytic'; end
    hf=figure; semilogy(sweepX.values, sweepX.(an).Psh_guided*1e6,'Linewidth',2); grid on;
    hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
    semilogy(sweepX.values, sweepX.(an).Psh_ScattAndGuided*1e6,'--','Linewidth',2);  
%     hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
%     semilogy(sweepX.values, sweepX.(an).Psh_ScattONLY*1e6,'-.','Linewidth',3);  
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('P_S_H_G (avg) [uW] ');
    legend([legnames]);  title({['P_S_H (',type,')'];replace(sweepX.fieldname,'_',' ')});
    % xlim([0,5000]);   
end
%--------------- Plot: EFFICIENCY
function makeplot_Efficiencies(legnames,sweepX,an)
    if strcmp(an,'n') type='numeric'; else type='analytic'; end
    hf=figure; loglog(sweepX.values, sweepX.(an).Effic_guided,'Linewidth',2); grid on;
    hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
    loglog(sweepX.values, sweepX.(an).Effic_ScatterAndGuided,'--','Linewidth',2);  
%     hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
%     loglog(sweepX.values, sweepX.(an).Effic_ScattONLY,'-.','Linewidth',3);  
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
    legend([legnames]);   title({['Efficiencies (',type,')'];replace(sweepX.fieldname,'_',' ')});
    % xlim([0,5000]); 
end
%--------------- Plot: (Scatt+Guided) / Guided
function makeplot_Scatt2Guided(legnames,sweepX,an)
    if strcmp(an,'n') type='numeric'; else type='analytic'; end
    
    hf=figure; semilogy(sweepX.values, cummax(sweepX.(an).Psh_ScattAndGuided)./cummax(sweepX.(an).Psh_guided) ,'Linewidth',2); grid on;
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('(Scatt+Guided) / Guided]');
    legend([legnames]);  title({['(Scatt+Guided) / Guided', '(',type,')'];replace(sweepX.fieldname,'_',' ')});
    xlim([min(sweepX.values),max(sweepX.values)]); 
end
%--------------- Plot: Coherence Length
function  makeplot_Coherence_Length(legnames,sweepX,an)
    if strcmp(an,'n') type='numeric'; else type='analytic'; end
    figure; plot(sweepX.values, sweepX.(an).Lcoh_um,'Linewidth',2); grid on;
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('Coherence Length [um]');
    legend(legnames);  title({['Coherence Length', '(',type,')'];replace(sweepX.fieldname,'_',' ')});    
end
 

%--------------- ScattONLY - Analytic v Numeric 
function makeplot_ScattONLY_Power_CompareAnalyticNumeric(legnames, sweepX)
    hf=figure; semilogy(sweepX.values, sweepX.a.Psh_ScattONLY*1e6,'Linewidth',2); grid on;
    hold on;  ax = gca; ax.ColorOrderIndex = 1;             %reset color ordering
    semilogy(sweepX.values, sweepX.n.Psh_ScattONLY*1e6,'--','Linewidth',2); ;  
    xlabel(replace(sweepX.fieldname,'_',' ')); ylabel('P_S_H_G (avg) [uW] - SCATT ONLY');
    legend([legnames]);  title({['P_S_H ScattONLY: Analytic(solid) v. Numeric(dashed)'];replace(sweepX.fieldname,'_',' ')});
    % xlim([0,5000]);   
end

function modan = apply_DF_to_results(modan)
    modan.Psh_guided            = modan.Psh_guided * modan.DF;
    modan.Psh_ScattONLY         = modan.Psh_ScattONLY * modan.DF;
    modan.Psh_ScatterAndGuided  = modan.Psh_ScatterAndGuided * modan.DF;
end


%% Calc NUMERIC
function modn = calcNumeric(c,IN)
    modn = IN;
    modn = calc_prefactor(c,modn);         %modn  holds data for *NUMERIC* integration result   
    % modn = calc_numeric_zerolinewidth(c,modn,delta_lam_nm,c.do_plot);  %CALL NUMERIC INTEGRATION over Z
    modn = calc_Numeric_WithLineWidth(c,modn);
    if modn.DF ~=1  modn = apply_DF_to_results(modn); end
end

%% OVERLAP:  Pump Linewidth X SHG Bandwidth  -- NUMERIC
function modn = calc_Numeric_WithLineWidth(c,modn)
    %Coherence Length (of PumpLinewidth)
    modn.Lcoh_um = (2 * 1.39156) /(2*pi) * (modn.lam1/modn.linewidth_nm) / modn.dispersion_slope_diff * 1e-3;
    
    if  modn.Lcoh_um > 2*modn.L_waveguide_um  %SKIP OVERLAP CALCULATION (L_coh > L_wg)       
        disp(['numeric : SKIP OVERLAP CALC:  L_coh=',num2str(round(modn.Lcoh_um,0)),' [um]  >  L_waveguide=',num2str(modn.L_waveguide_um),' [um]']);
        delta_lam_nm = modn.Pump_Detune_nm;
        modn = calc_numeric_zerolinewidth(c,modn,delta_lam_nm);         %CALL NUMERIC INTEGRATION over Z
        modn.Overlap_Applied           = 'No';        
        modn.Ppump_remain              = modn.Ppump_remain;        %no overlap (zero linewidth)
        modn.Psh_guided                = modn.Psh_guided;          %no overlap (zero linewidth)      
        modn.Psh_ScattONLY             = modn.Psh_ScattONLY;       %no overlap (zero linewidth)  
        modn.Psh_ScatterAndGuided      = modn.Psh_ScatterAndGuided ;%no overlap (zero linewidth)
        
    else     %DO OVERLAP CALCULATION (L_coh > L_wg)  
        
        %Pump Linewidth 
        N_steps_integration = c.N_steps_overlap_integration;
        HWHM = modn.linewidth_nm/2;
        deltas_lam_nm = linspace(-2*HWHM-modn.Pump_Detune_nm,2*HWHM+modn.Pump_Detune_nm, N_steps_integration )'; %[nm]
        Ypump = exp(-0.5*((deltas_lam_nm - modn.Pump_Detune_nm).^2)/(HWHM)^2 );
        
        %DO OVERLAP INTEGRATION:   
        for klam = 1:length(deltas_lam_nm) 
            modn = calc_numeric_zerolinewidth(c,modn,deltas_lam_nm(klam) );   %CALL NUMERIC INTEGRATION over Z
            Ppump_remain(klam,1)            = modn.Ppump_remain;
            Psh_guided(klam,1)              = modn.Psh_guided;
            Psh_ScattONLY(klam,1)           = modn.Psh_ScattONLY;
            Psh_ScatterAndGuided(klam,1)    = modn.Psh_ScatterAndGuided;
        end
        modn.Overlap_Applied           = 'Yes';
        modn.Ppump_remain              = sum(Ypump.*Ppump_remain)/sum(Ypump);         % **** numeric overlap integration here** 
        modn.Psh_guided                = sum(Ypump.*Psh_guided)/sum(Ypump);           % **** numeric overlap integration here** 
        modn.Psh_ScattONLY             = sum(Ypump.*Psh_ScattONLY)/sum(Ypump);        % **** numeric overlap integration here** 
        modn.Psh_ScatterAndGuided      = sum(Ypump.*Psh_ScatterAndGuided)/sum(Ypump); % **** numeric overlap integration here** 
        
        
        %*** plot: Linewidth(Pump) v. Bandwidth(PhaseMatching)
        if c.do_plots==1 %|  (modn.Efficiency_ScatterAndGuided <  modn.Efficiency)
            figure; subplot(2,1,1); plot(deltas_lam_nm,Ypump); 
            hold on; yyaxis right; plot(deltas_lam_nm,[Psh_guided,Psh_ScattONLY,Psh_ScatterAndGuided]);    ylabel('Abs. Value of POWER');       
            legend('Pump','guided',...
                   'ScattONLY',...
                   'guided+Scatt');
            title({'OVERLAP: PhaseMatch Bandwidth v. Pump Linewidth  (*Numeric*)';...
                ['Overlap Power(guided+scatt) = ',num2str(modn.Psh_ScatterAndGuided)]});
            
            subplot(2,1,2);
            plot(deltas_lam_nm,Ypump);  hold on; yyaxis right; 
            plot(deltas_lam_nm,[Psh_guided/max(Psh_guided),Psh_ScattONLY/(max(Psh_ScattONLY)),...
                        Psh_ScatterAndGuided/max(Psh_ScatterAndGuided)]);   
            ylabel('Rescaled each Signal to 1'); 
            xlabel('\Delta WL [nm]');
        end
        
        
        
    end %endif:  do overlap integral
end

%%  Calc Numeric:  ZERO LINEWIDTH  (With Detuning,  With Loss)
function  modn = calc_numeric_zerolinewidth(c,modn,delta_lam_linewidth_nm)

       %run NUMERIC Integation  ..................................
        
        %parameters define:
        delta_lam = modn.Pump_Detune_nm + delta_lam_linewidth_nm;     %this is delta from ideal PhaseMatch WL.  it is sum of:
                                                                 % Detuning (Center of Pump Line) +  Offset from Linewidth (modn.linewidth_nm/2) * 2;
        delta_n           = delta_lam * modn.dispersion_slope_diff;    %unitless
        delta_k_cm        = 4*pi/(modn.lam1*1e-7) * delta_n;           %[cm^-1]
        modn.alpha0_cm    = modn.Loss_pump_dBcm / 4.3429;         % pump:          alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
        modn.alpha3_cm    = modn.Loss_SH_dBcm / 4.3429;           % SH:            alpha[cm-1] =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
        modn.alpha3abs_cm = modn.Loss_abs_SH_dBcm / 4.3429;       % SH: ABSorption alpha[cm-1] =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
        dk_um             = delta_k_cm     * (1e-4);
        a0_um             = modn.alpha0_cm  * (1e-4); 
        a3_um             = modn.alpha3_cm  * (1e-4); 
        g                 = modn.g;
        a3a_um            = modn.alpha3abs_cm  * (1e-4); 

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
        dP3_loss = 0; dP3_gain=0;
        modtmp=modn;
        for kz = 1:length(zvec)-1   
            
            %track change in P3guided
            dP3_loss(kz+1,1) = (-0.5*a3_um*A3(kz,1)*dz_um)*conj(-0.5*a3_um*A3(kz,1)*dz_um);
            dP3_gain(kz+1,1) = (i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz))*dz_um)...
                            *conj(i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz))*dz_um);
                        
            %NO ABSORPTION -- SH Loss = ALL SCATTER
            %diffeq_1:    dA3/dz = -0.5a3*A3  + i*g*A1^2*exp(-i*dk*dz)
%             dA3 = ( -0.5*a3_um*A3(kz,1)          +  i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz)) ) * dz_um;
%             A3(kz+1,1) = A3(kz,1) + dA3;
            
            %WITH ABSORPTION -- SH Loss = ABS + SCATTER
            %diffeq_1:    dA3/dz = -0.5a3*A3  + i*g*A1^2*exp(-i*dk*dz)
            dA3 = ( -0.5*(a3_um+a3a_um)*A3(kz,1) +  i * g * abs(A1(kz,1))^2 *exp(-i*dk_um*zvec(kz)) ) * dz_um;
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
            
                       
            %for comparison only:   analytic result:  def_Integral(0...L)
            L                       = zvec(kz);
            modtmp.L_waveguide_um   = L;
            modtmp = calc_prefactor(c,modtmp); 
            modtmp = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,modtmp);        %
            modtmp = PhaseMismatch_and_Loss_Derate(c, modtmp, modtmp.Pump_Detune_nm);   % analytic solution            
            Eff_zlwnl(kz+1,1)               = modtmp.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss;
            Eff_guided(kz+1,1)              = modtmp.DeRate_Wei * Eff_zlwnl(kz+1,1);              
            Eff_scatt(kz+1,1)               = modtmp.DeRate_Scatt * Eff_zlwnl(kz+1,1);               
            P3z_analytic(kz+1,1)            = modtmp.Pp * Eff_guided(kz+1,1);           %[W]
            PscattONLY_analytic(kz+1,1)     = modtmp.Pp * Eff_scatt(kz+1,1);            %[W]
            P3plusScatt_analytic(kz+1,1)    = modtmp.Pp * Eff_scatt(kz+1,1);            %[W]
 
        end %Integrate_over_Z
        
%         figure; plot(zvec,[dP3_gain,dP3_loss]) 

        % Amplitude --> Power
        P1z                 = A1.*conj(A1);                                   %[W]
        P3z                 = A3.*conj(A3);                                   %[W]
        P3scatt_z           = A3scatt.*conj(A3scatt);                         %[W]
        %Value at z=L
        modn.Ppump_remain           = P1z(end);
        modn.Psh_guided             = P3z(end);
        
        %amplitude method
%         modn.Psh_ScattONLY          = P3scatt_z(end); %original method
%         modn.Psh_ScatterAndGuided 	= P3z(end)+P3scatt_z(end);

        %Wei's  method
        modn.Psh_ScattONLY          = P3b_scatt_z(end); %wei's method        
        modn.Psh_ScatterAndGuided 	= P3b_scatt_z(end)+P3z(end);
        
   
        % COMPARE:  NUMERICAL v. ANALYTICAL
        if c.do_compare_to_analytic==1
            if delta_lam == 0
                disp([P3z_analytic(end), modn.Psh_guided]);
            end                    
            figure; plotyy(zvec, P1z, zvec, [P3z, P3scatt_z, P3b_scatt_z] ); ylabel('power [W]'); xlabel('z [um]');
            legend('P1','P3','P3scatt','P3bscatt(Wei)');    title({'Numeric Results (only)',['delta WL = ',num2str(delta_lam),'[nm]']});
            figure; semilogy(zvec, [P1z, P3z, P3scatt_z, P3b_scatt_z] ); ylabel('power [W]'); xlabel('z [um]');
            legend('P1','P3','P3scatt','P3bscatt(Wei)');    title({'Numeric Results (only)',['delta WL = ',num2str(delta_lam),'[nm]']});
            
            figure;semilogy(zvec, [P3z, P3z_analytic, P3scatt_z, PscattONLY_analytic]);  ylabel('power [W]'); xlabel('z [um]');
            legend('guided: numeric','guided: analytic (no detune)','scatt(Only): numeric','scatt(Only): analytic (no detune)');  
            title({'Compare:  Numeric v. Analytic ',['delta WL = ',num2str(delta_lam),'[nm]']});
        end        

end  %NUMERIC CALC






%%  Calc ANALYTIC
function moda = calcAnalytic(c,IN)
    moda = IN;
    moda = calc_prefactor(c,moda); 
    moda = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,moda);
    moda = DeRate_NonZeroLinewidth_analytic(c,moda);    
    if moda.DF ~=1  moda = apply_DF_to_results(moda); end
end


%%  Prefactor Calc
function moda = calc_prefactor(c,moda) 
    moda.lam2            = moda.lam1 / 2;
%     moda.w_1             = 2*pi* (c.c0/moda.n1) / (  moda.lam1 * 1e-9);    %[Hz]
%     moda.w_2             = 2*pi* (c.c0/moda.n2) / (  moda.lam2 * 1e-9);    %[Hz]    
    moda.w_1             = 2*pi* (c.c0) / (  moda.lam1 * 1e-9);    %[Hz]
    moda.w_2             = 2*pi* (c.c0) / (  moda.lam2 * 1e-9);    %[Hz]  
    moda.prefactor       = 2*moda.w_2^2*(moda.d_xx*1e-12)^2 / (c.eps0*moda.n1^2*moda.n2*c.c0^3);
    moda.A_waveguide_um2 = moda.x_waveguide_um * moda.y_waveguide_um;     %[um^2]
    moda.g               = sqrt( moda.overlap_frac^2 / moda.A_waveguide_um2 * moda.prefactor );           
end

%%
function moda = Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss(c,moda)
    moda.A_waveguide_um2  = moda.x_waveguide_um * moda.y_waveguide_um;  %[um^2]
    moda.GAMMA_pctPerWatt = moda.prefactor * moda.overlap_frac^2 ...
                            * (moda.L_waveguide_um^2 / moda.A_waveguide_um2);  
    moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss = moda.GAMMA_pctPerWatt * moda.Pp;
end

%%
function moda = PhaseMismatch_and_Loss_Derate(c, moda, delta_lam_nm)
    delta_lam   = delta_lam_nm;                             %(moda.linewidth_nm/2) * 2;
    delta_n     = delta_lam * moda.dispersion_slope_diff;    %unitless
    delta_k_cm  = 4*pi/(moda.lam1*1e-7) * delta_n;           %[cm^-1]
    moda.alpha0_cm = moda.Loss_pump_dBcm / 4.3429;           % pump: alpha[cm-1] =  [dB/cm] * (ln(10)/10) = [dB/cm]/4.34
    moda.alpha3_cm = moda.Loss_SH_dBcm / 4.3429;             % SH: alpha[cm-1]   =  [dB/cm] * (ln(10)/10  = [dB/cm]/4.34
    
    dkL  = delta_k_cm   * (moda.L_waveguide_um * 1e-4);
    a3L   = moda.alpha3_cm * (moda.L_waveguide_um * 1e-4); 
    a0L   = moda.alpha0_cm * (moda.L_waveguide_um * 1e-4); 
    daL   = a3L/2 - a0L;    
    
    %*** GUIDED LIGHT (Only) ****
    if dkL==0   
        moda.DeRateAppx0     = (exp(-a3L));
        moda.DeRateAppx      = (exp(-a3L));
        moda.DeRateAppx2     = (exp(-a3L));
        moda.DeRate_BF       = (exp(-a3L)) *2* ( cosh(a3L) - 1 ) / ( (a3L)^2 ) ;
        if daL==0
            moda.DeRate_Wei  = (exp(-a3L))*sinc(dkL/2/pi)^2;
        else
            moda.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
                                          / ( (dkL)^2 + (daL)^2 ) ;
        end
        %following 4 variables are here just to put them into struct
        moda.DeRate_Scatt            = []; 
        moda.DeRate_Scatt_indefInteg = [];
        moda.DeRate_Wei_IdealSHLoss  = [];
        moda.DeRate_Scatt_alt1       = [];
        if a0L==0 
            moda.DeRate_noPMreq  = 1;       
        else
            %original formula (brent's "incorrect" version):
            moda.DeRate_noPMreq  = (exp(-a0L)-1)^2 / ( (a0L)^2 ) ;       
        end      
        
    else %dkL>0
        % crude estimate based on SINC(x) only ....
        moda.DeRateAppx0 = (exp(-a3L)) *1* ( sinc(dkL/2) )^2 ;
        moda.DeRateAppx  = (exp(-a3L)) *1* ( sin(dkL/2) )^2 / (dkL/2)^2  ;
        moda.DeRateAppx2 = (exp(-a3L)) *1* ( sin(dkL/2) /  (dkL/2)  )^2  ;   
        
        %original formula (brent's "incorrect" version):
        moda.DeRate_BF     = (exp(-a3L)) *2* ( cosh(a3L) - cos(dkL) ) / ( (dkL)^2 + (a3L)^2 ) ;

        %formula from Wei:
        moda.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
                                          / ( (dkL)^2 + (daL)^2 ) ;   
        
        moda.DeRate_noPMreq = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*1 ) ...
                                          / ( (daL)^2 ) ;                                          
    end
      
    %*** SCATTERED LIGHT (Only) ***
    moda.DeRate_Scatt    = ((0.5*a3L)^2) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
                            % this results from using *DEFINITE*
                            % integration in the derivation
                            
    moda.DeRate_Scatt_indefInteg    ...
                            = ((0.5*a3L)^2) * ( exp(-2*a0L) ) ...
                            / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
                            % note: no sinc() function here -- (phasematch relaxed?)
                            % this results from using *INdefinite*
                            % integration in the derivation                  
                        
    moda.DeRate_Wei_IdealSHLoss  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / ( (dkL)^2 + (a0L)^2 );
                            %this sets alpha_SH Loss = 0 
                            % .... but still has phase matching requirement
            
    moda.DeRate_Scatt_alt1  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
                            / ( (dkL)^2 + (a0L)^2 )   -  moda.DeRate_Wei;
                            % "Wei_Ideal  - Wei"
                        
  
end

%% OVERLAP: Pump Linewidth X SHG Bandwidth
function moda = DeRate_NonZeroLinewidth_analytic(c,moda)

    %Coherence Length (of PumpLinewidth)
    moda.Lcoh_um = (2 * 1.39156) /(2*pi) * (moda.lam1/moda.linewidth_nm) / moda.dispersion_slope_diff * 1e-3;
    
    
    if  moda.Lcoh_um > 2*moda.L_waveguide_um  %SKIP OVERLAP CALCULATION (L_coh > L_wg)        
        disp(['analytic: SKIP OVERLAP CALC:  L_coh=',num2str(round(moda.Lcoh_um,0)),' [um]  >  L_waveguide=',num2str(moda.L_waveguide_um),' [um]']);
        delta_lam_nm = moda.Pump_Detune_nm;
        moda = PhaseMismatch_and_Loss_Derate(c,moda, delta_lam_nm);
        
        moda.DeRate_Overlap_guided       = moda.DeRate_Wei;                        % Spectral Overlap = 1 
        moda.DeRate_Overlap_ScattONLY    = moda.DeRate_Scatt;                      % Spectral Overlap = 1         
        moda.Efficiency_guided           = moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss * moda.DeRate_Overlap_guided;    
        moda.Efficiency_ScattONLY        = moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss * moda.DeRate_Overlap_ScattONLY; 
        moda.Efficiency_ScatterAndGuided = moda.Efficiency_guided  + moda.Efficiency_ScattONLY;                
        moda.Psh_guided                  = moda.Pp * moda.Efficiency_guided;
        moda.Psh_ScattONLY               = moda.Pp * moda.Efficiency_ScattONLY; 
        moda.Psh_ScatterAndGuided        = moda.Pp * moda.Efficiency_ScatterAndGuided;

        
    else  %(moda.Lcoh_um < 2*moda.L_waveguide_um  -->  DO OVERLAP CALCULATION

        %Pump Linewidth 
        HWHM = moda.linewidth_nm/2;
        N_steps_integration = c.N_steps_overlap_integration;
        deltas_lam_nm = linspace(-2*HWHM-+moda.Pump_Detune_nm,2*HWHM+moda.Pump_Detune_nm, N_steps_integration )'; %[nm]
        Ypump = exp(-0.5*((deltas_lam_nm - moda.Pump_Detune_nm).^2)/(HWHM)^2 );

        %Bandwidth of Phase Matching -- generate phasematching de-rate curve (over domain of 2x linewidth) 
        for k = 1:length(deltas_lam_nm)  
            moda = PhaseMismatch_and_Loss_Derate(c,moda,deltas_lam_nm(k));
                DeRate_BF(k,1)      = moda.DeRate_BF;
                DeRateAppx(k,1)     = moda.DeRateAppx;
                DeRate_noPMreq(k,1) = moda.DeRate_noPMreq;
                DeRate_Wei(k,1)     = moda.DeRate_Wei ;
                DeRate_ScattONLY(k,1)   = moda.DeRate_Scatt;               
        end
        DERATE_TO_USE = DeRate_Wei;  % <----  !!! select which expression to use  !!!!
    %     DERATE_TO_USE = DeRate_noPMreq;  % <----  !!! select which expression to use  !!!!

        %OUTPUT:  Overlap:  Linewidth(pump) x Bandwidth(Phasematch)
        moda.DeRate_Overlap_guided       = sum(Ypump.*DeRate_Wei)/sum(Ypump);                 % **** numerical integration here** 
        moda.DeRate_Overlap_ScattONLY    = sum(Ypump.*DeRate_ScattONLY)/sum(Ypump);           % **** numerical integration here** 
        
        modn.Overlap_Applied            = 'Yes';
        moda.Efficiency_guided           = moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss * moda.DeRate_Overlap_guided;    
        moda.Efficiency_ScattONLY        = moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss * moda.DeRate_Overlap_ScattONLY; 
        moda.Efficiency_ScatterAndGuided = moda.Efficiency_guided + moda.Efficiency_ScattONLY;
        moda.Psh_guided                  = moda.Pp * moda.Efficiency_guided;
        moda.Psh_ScattONLY               = moda.Pp * moda.Efficiency_ScattONLY; 
        moda.Psh_ScatterAndGuided       = moda.Pp * moda.Efficiency_ScatterAndGuided;
        
        
        if isnan(moda.Efficiency_guided)
            warning('moda.Efficiency is NAN ')
            disp(moda);
        end

        %*** plot: Linewidth(Pump) v. Bandwidth(PhaseMatching)
        if c.do_plots==1 %|  (moda.Efficiency_ScatterAndGuided <  moda.Efficiency)
            figure; subplot(2,1,1);
            plot(deltas_lam_nm,Ypump);
            hold on; yyaxis right; plot(deltas_lam_nm,[DeRate_Wei,  DeRate_ScattONLY]);  ylabel('Abs. Value of DeRate');        
            legend('Pump Linewidth',...
                   'DeRate: PhaseMatching (Wei)',...
                   'DeRate: Scatt Only');
            title({'OVERLAP:  PhaseMatch Bandwidth v. Pump Linewidth  (*Analytic*)';...
                ['Overlapped Derate = ',num2str(moda.DeRate_Overlap_guided)]});
      
            subplot(2,1,2);
            plot(deltas_lam_nm,Ypump);  hold on; yyaxis right; 
            plot(deltas_lam_nm,[DeRate_Wei/max(DeRate_Wei),  DeRate_ScattONLY/max(DeRate_ScattONLY)]);   
            xlabel('\Delta WL [nm]');  ylabel('Rescaled');
        end
               
    end %endif:  do overlap integral

    
end  %end:  DeRate_NonZeroLinewidth_analytic







%% get INPUT 
function INs = get_INPUT(dirs,inputfile,sheetname)
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
            INs{kc}.Loss_abs_SH_dBcm = tmp.num(15,kc); 
            INs{kc}.Pump_Detune_nm   = tmp.num(16,kc); 
        end
        Ncases = length(INs);
        
    else
        error('Excel Input was not obtained... Returning Default INPUTs ')
    end
    
end


%% Generate Output Matrix
function M_out = ResultsMatrix_AddColumn(MR, moda, modn)
   newres      = zeros(5,1);
   newres(1,1) = moda.Pp;
   newres(2,1) = moda.Lcoh_um;
   newres(3,1) = moda.prefactor;
   newres(4,1) = moda.GAMMA_pctPerWatt; 
   newres(5,1) = moda.Efficiency_ZeroLinewidth_and_NoDetune_and_NoLoss;
   newres(6,1) = moda.DeRate_Overlap_guided;
   newres(7,1) = moda.Psh_guided;
   newres(8,1) = moda.Psh_ScattONLY;
   newres(9,1) = moda.Psh_ScatterAndGuided;
   newres(10,1) = -99;  %marker for Numeric results   
   newres(11,1) = -99;  %marker for Numeric results
   %Numeric Results
   newres(12,1) = modn.Psh_guided;
   newres(13,1) = modn.Psh_ScattONLY;
   newres(14,1) = modn.Psh_ScatterAndGuided;
   newres(15,1) = modn.Ppump_remain;
   if isfield(modn,'Overlap_Applied') 
       if strcmp(modn.Overlap_Applied,'Yes')
           newres(16,1) = 1; 
       else
           newres(16,1) = 0;
       end
   else
       newres(16,1) = 0;
   end
   
   
    
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



%%  TEST DERATE FUNCTION(S)

function test_DeRate_Functions(c,moda)

    delta_lam_nm_list = linspace(0,0.5,101)';
    a0_list = linspace(0,20,11)';
    a3_list = linspace(0,1000,51)';
    k1 = 1;
    moda.alpha0_cm = a0_list(1);
    clear D D2
    close all;
    moda.L_waveguide_um = 200;
    for k1=1:length(delta_lam_nm_list)
    %     for k0=1:length(a0_list)
    %         moda.Loss_pump_dBcm = a0_list(k0);
            for k3=1:length(a3_list) 
               moda.Loss_SH_dBcm = a3_list(k3)    
               moda          = PhaseMismatch_and_Loss_Derate(c, moda, delta_lam_nm_list(k1));           
               D(k3,k1)     = moda.DeRate_Wei;               %GUIDED
               D2(k3,k1)    = moda.DeRate_Scatt;             %scatter (only)
               D20(k3,k1)   = moda.DeRate_Scatt_indefInteg;  %scatter (only)
               D3a(k3,k1)   = moda.DeRate_Scatt_alt1;        %scatter (only)
               D3b(k3,k1)   = moda.DeRate_Wei_IdealSHLoss;   %GUIDED  

               Dsum(k3,k1)  = moda.DeRate_Wei +  moda.DeRate_Scatt;
            end
    %     end
    end

     figure; imagesc(delta_lam_nm_list,a3_list,log10(D(:,:)));  ylabel('Loss(SH) [dB/cm]'); xlabel('WL mismatch [nm]');  
     title({['DWei (guided!);   Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    colorbar; caxis([-8,0]);
     figure; imagesc(delta_lam_nm_list,a3_list,log10(D2(:,:)));  ylabel('Loss(SH) [dB/cm]'); xlabel('WL mismatch [nm]'); 
     title({['Dscatt(orig,def); Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    colorbar; caxis([-8,0]);
     figure; imagesc(delta_lam_nm_list,a3_list,log10(D20(:,:)));  ylabel('Loss(SH) [dB/cm]'); xlabel('WL mismatch [nm]'); 
     title({['DeRate.Scatt.IndefInteg(orig,indef); Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    colorbar; caxis([-8,0]);
     figure; imagesc(delta_lam_nm_list,a3_list,log10(Dsum(:,:)));  ylabel('Loss(SH) [dB/cm]'); xlabel('WL mismatch [nm]'); 
     title({['DWei + Dscatt(def.integral);   Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    colorbar; caxis([-8,0]);


    figure; plot(delta_lam_nm_list,log10(D(:,:)));   xlabel('WL mismatch [nm]');  ylabel('Derate(Wei)');  ylim([-10,0]);
      title({['DWei;                      Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    figure; plot(delta_lam_nm_list,log10(D2(:,:)));   xlabel('WL mismatch [nm]');  ylabel('Dscatt(orig,def.integral)'); ylim([-10,0]);
      title({['Dscatt(orig,def.integral); Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    figure; plot(delta_lam_nm_list,log10(D20(:,:)));   xlabel('WL mismatch [nm]');  ylabel('Dscatt(orig,Indef.integral)'); ylim([-10,0]);
      title({['Dscatt0( INDEF.integral); Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});
    figure; plot(delta_lam_nm_list,log10(Dsum(:,:)));   xlabel('WL mismatch [nm]');  ylabel('Derate  '); ylim([-10,0]);
      title({['DWei + Dscatt(def.integral);   Loss (at Pump): ',num2str(moda.Loss_pump_dBcm),'dB/cm'];['Length: ',num2str(moda.L_waveguide_um),'um']});

end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD -- NOT USED 
%% Derate (Simple)
% function TestDerateSimple()  %%mimic Excel output
%     losses0 = linspace(2,20,10);
%     losses3 = linspace(2,20,10);
%     for k0 = 1:length(losses0)
%         for k1 = 1:length(losses3)
%     end
%         
% end

% function out = DerateSimple(alpha0_cm, alpha0_cm, delta_k_cm, L_waveguide_um)
%     dkL  = delta_k_cm    * (L_waveguide_um * 1e-4);
%     a3L  = alpha3_cm     * (L_waveguide_um * 1e-4); 
%     a0L  = moda.alpha0_cm * (L_waveguide_um * 1e-4); 
%     daL  = a3L/2 - a0L; 
%     
%     %use this to write code that i can verify matches Excel output
%     
%         %formula from Wei:
%         out.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
%                                       / ( (dkL)^2 + (daL)^2 ) ;        
%         %SCATTER SUM
%         out.DeRate_Scatt = ((0.5*a3L)^2) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%                             / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  );                        
%         out.DeRate_Scatt_indefInteg    = ((0.5*a3L)^2) * ( exp(-2*a0L) ) ...
%                             / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
%         out.DeRate_Wei_IdealSHLoss  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%                             / ( (dkL)^2 + (a0L)^2 );
%                                                          
%         out.DeRate_Scatt_alt1  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%             / ( (dkL)^2 + (a0L)^2 )   -  out.DeRate_Wei;
%     
%     
%     
% end




 %from  DERATE Function ...       
%     else %dkL>0
%         % crude estimate based on SINC(x) only ....
%         moda.DeRateAppx0 = (exp(-a3L)) *1* ( sinc(dkL/2) )^2 ;
%         moda.DeRateAppx  = (exp(-a3L)) *1* ( sin(dkL/2) )^2 / (dkL/2)^2  ;
%         moda.DeRateAppx2 = (exp(-a3L)) *1* ( sin(dkL/2) /  (dkL/2)  )^2  ;
%        
%         %original formula (brent's "incorrect" version):
%         moda.DeRate_BF     = (exp(-a3L)) *2* ( cosh(a3L) - cos(dkL) ) / ( (dkL)^2 + (a3L)^2 ) ;
% 
%         %formula from Wei:
%         moda.DeRate_Wei  = (exp(-a3L)) * ( (exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2 ) ...
%                                       / ( (dkL)^2 + (daL)^2 ) ;
%         
%         %SCATTER SUM
%         moda.DeRate_Scatt = ((0.5*a3L)^2) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%                             / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  );
%                             % this results from using *DEFINITE*
%                             % integration in the derivation
%                         
%         moda.DeRate_Scatt_indefInteg   ...
%                             = ((0.5*a3L)^2) * ( exp(-2*a0L) ) ...
%                             / (  ( (dkL)^2 + (daL*a0L) )^2   + ( (daL*dkL)-(a0L*dkL) )^2  ); 
%                             % note: no sinc() function here -- (phasematch relaxed?)
%                             % this results from using *INdefinite*
%                             % integration in the derivation        
% 
%         moda.DeRate_Wei_IdealSHLoss  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%                             / ( (dkL)^2 + (a0L)^2 );
%                             %this sets alpha_SH Loss = 0 
%                             % .... but still has phase matching requirement
%                             
%                                                          
%         moda.DeRate_Scatt_alt1  = (exp(0)) * ( (exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2 ) ...
%                             / ( (dkL)^2 + (a0L)^2 )   -  moda.DeRate_Wei;
%                             
%                         
%     end 
%     
%     
%         if isnan(moda.DeRate_Scatt)
%             disp(moda);
%         end
%         if isnan(moda.DeRate_Wei)
%             disp(moda);
%         end
%         
        
    
%     if dkL>0
%         moda.DeRate_BF=  2 * ( cosh(aL) - cos(dkL) ) * exp(-aL) / ( (dkL)^2 + (aL)^2 ) ;
%         moda.DeRateAppx = 2 *( sinc(dkL) )^2 * (exp(-aL));
%     else
%         if aL>0
%             moda.DeRate_BF=  exp(-aL);
%             moda.DeRateAppx =  exp(-aL);
%         else
%             moda.DeRate_BF= 1;
%             moda.DeRateAppx = 1;
%         end
%     end
    





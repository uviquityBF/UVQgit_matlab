function process_FLIR_CAL_Data2_manyTexp()

    clear;
    close all;
    dirs.npy = '/Volumes/GoogleDrive/Shared drives/Engineering/Analysis (Brent3)/MATLAB/npy-matlab';
%     dirs.input = '/Volumes/GoogleDrive/Shared drives/Engineering/Optical Hazard - Eye Safety/FLIR_camera_setup_and_calibration/2020_06_07_flir_calibration_data';
    dirs.input = '/Volumes/GoogleDrive/Shared drives/Engineering/Optical Hazard - Eye Safety/FLIR_camera_setup_and_calibration/';
    dirs.output = [dirs.input,'output'];
 
    %load images
    cd(dirs.input);
    [filenames,path,dummy] = uigetfile('*.mat','select file','MultiSelect','on');
    
    Nfiles = length(filenames);
    
   %sort images by gain
    for kf = 1:Nfiles
        s0 = split(filenames{kf},{'gain_','Gain_','Gain'});
        s1 = split(s0{2},{'_time_','_Time','_Time_'});
        gainval(kf,1) = str2num(s1{1});
    end
   [sorted_gain, indc_sort] = sort(gainval);
   gainval = sorted_gain;
   
   
   %get RadiativePower
   default_power = 0.0003979; %[W]
   default_power = 0.00028; %[W]
   default_power = 0.000066; %[W]
   default_power = 0.000075; %[W]
   ans = inputdlg('What was the measured (time-average) radiant power?  [W]','True Power',1,{num2str(default_power)});
   Power_W = str2num(ans{1});
   
   
   % ------- TIF FILES -------
   tif_sat_thresh_DN = 255; 
   if strcmp(filenames{1}(end-2:end),'tif')  
        for kf = 1:Nfiles
            IMGs{kf}  = imread([path,filenames{indc_sort(kf)}]);
            s0 = split(filenames{indc_sort(kf)},{'gain_','Gain'});
            s1 = split(s0{2},'_time_');
            gain_check(kf,1) = str2num(s1{1});
            s2 = split(s1{2},'_');
            Texp_us(kf,1)   = str2num(s2{1});  
            SatMats{kf}     = (IMGs{kf} >=tif_sat_thresh_DN);  %pixels saturated           
            isel_saturated  = find(SatMats{kf});
            NumSatPix(kf,1) = length(isel_saturated(:));
        end
       
   else %  ------- MAT FILES -------  ... from pydevcam collection
       %read in file data in Sorted Order  
        for kf = 1:Nfiles
            IN{kf} = load([path,filenames{indc_sort(kf)}]);
            s0 = split(filenames{indc_sort(kf)},{'gain_','Gain_','Gain'});
            s1 = split(s0{2},{'_time_','_Time','_Time_'});
            gain_check(kf,1) = str2num(s1{1});
            s2 = split(s1{2},{'_','us_'}); if length(s2{1})==0 s2{1}=s2{2}; end
            Texp_us(kf,1) = str2num(s2{1});              
            NumSatPix(kf,1) = sum(IN{kf}.imagedict{1}.saturation(:));
            SatMats{kf} = IN{kf}.imagedict{1}.saturation;
            
            %get pixel map with no saturated pixels
            IMGs{kf} = IN{kf}.imagedict{1}.image;
        end

   end

   
   IMG0 = IMGs{1};    %UNSATURATED IMAGE  (for estimation of missing power in saturated images)    

   %compute J/count for each FIle  & Flag Saturatin
   for kf = 1:Nfiles        
       
        %use *UNSATURATED IMAGE* to make estaimte of how much power is
        %missing from any set of saturated pixels ....
        %remove saturated pixels from Sum
        satPix = find(SatMats{kf}(:)==1);
        frac_sum_in_satPix(kf,1) = sum( IMG0(satPix) ) / sum(IMG0(:));
        %--------end:  estimation of fraction of power/sum in saturated pixels
        
        
        %sum_ with NO Sat. Pixels
     	validPix = find(SatMats{kf}(:)==0);
        SumCount_NoSatPix(kf,1) = sum( IMGs{kf}(validPix) );
       	SumCount_AllPix(kf,1) = sum( IMGs{kf}(:) );

        %Counts per Joule: AllPix
      	Joules_per_Image_AllPix(kf,1)   = Power_W * (Texp_us(kf,1)*1e-6);
        Count_per_J_AllPix(kf,1) =  SumCount_AllPix(kf,1) / Joules_per_Image_AllPix(kf,1);

        %Counts per Joule: NoSatPix
    	Joules_per_Image_NoSatPix(kf,1) = Power_W * (Texp_us(kf,1)*1e-6)*(1-frac_sum_in_satPix(kf,1));
     	Count_per_J_NoSatPix(kf,1) =  SumCount_NoSatPix(kf,1) / Joules_per_Image_NoSatPix(kf,1);
       
        %SatFraction
        frac_pix_sat(kf,1) = length(satPix) / length(SatMats{kf}(:)); %fraction of pixels saturated
        
      	%generate figure of all images
      	nx = ceil( sqrt(Nfiles) );
    	ny = nx;
     	subplot(nx,ny,kf); imagesc(IMGs{kf}); colormap('bone');

   end
        
        

   frac_pix_sat_thresh = 0.1;   %thresh of allowed saturation in file to be used in fits 
   ans = inputdlg('In order for image to be used in calibration: \n Allow what fraction of Pixels to be saturated?','Frac Pixels Saturated',1,{num2str(frac_pix_sat_thresh)});
   frac_pix_sat_thresh= str2num(ans{1}); %thresh of allowed saturation in file to be used in fits 
    
   %2D Fits ------------------------------------------------------
    
    flag.flat = 2;

   iselfit = find(frac_pix_sat < frac_pix_sat_thresh  &   frac_sum_in_satPix < frac_pix_sat_thresh );
   if length(iselfit)==0  error('No data files meet sat_thresh requirement'); end
  
  
 %    figure; scatter3(xG,yT, SumCount_AllPix);
   figure; scatter3(gainval,Texp_us, SumCount_AllPix, 'o');
   hold on;  scatter3(gainval(iselfit),Texp_us(iselfit), SumCount_AllPix(iselfit),'x');
   zlabel('Sum Counts (All Pixels)');  ylabel('Gain'); xlabel('Texp [us]');
 
   %set up and run optimization
   xF = gainval(iselfit);  
   yF = Texp_us(iselfit); 
   zF = SumCount_AllPix(iselfit);
   z0 =  median(SumCount_AllPix(iselfit));
    %starting point
    ax = 1e8/20;    ay = 5;         %linear coefficients
    bx = 0;   by = 0;    bxy = 0;   %second order
    vars(1) = z0;
    vars(2) = ax;
    vars(3) = ay;
    vars(4) = bx;
    vars(5) = by;
    vars(6) = bxy;
    vars = double(vars);

    if flag.flat == 1
        function2optimize = @(vars) zDiff_function(vars,xF,yF,zF,@zFitpoly1_function);  %indicates which of the inputs to Diff function to optimize
    elseif flag.flat == 2
        function2optimize = @(vars) zDiff_function(vars,xF,yF,zF,@zFitpoly2_function);  %indicates which of the inputs to Diff function to optimize
    end
    optoptions = optimset('Display','iter','MaxIter',1000);
    [vars_bestfit, fval, exitflag, output] = fminsearch( function2optimize, vars, optoptions);
     
     %show fit:  Z
    nfit = 10;
    [x_show,y_show] = meshgrid(linspace(min(xF)-2,max(xF)+2,nfit),linspace(min(yF)-2,max(yF)+2,nfit)); 
    z_show = zFitpoly2_function(vars_bestfit,x_show,y_show);
    output.hf1 = figure; mesh(x_show,y_show,z_show); 
    hold on; plot3(xF,yF,zF,'o','Color','k','MarkerSize',2);
    legend('fit','data');   zlabel('Sum Counts (All Pixels)');  ylabel('Gain'); xlabel('Texp [us]');
    fitparams.twoD.params = vars_bestfit;
  
    % 1D Fits: v. Texp  (@ each gain)  ------------------------------------------------------
    gainvals_rounded = round(gainval,0);
    uGainVals = unique(gainvals_rounded)
    
    hf1 = figure; hf2=figure;
    for kG = 1:length( uGainVals )
        Gain_sel = uGainVals(kG);  legnames{kG} = [num2str(Gain_sel)];
        isel = find(Gain_sel == gainvals_rounded  ...
                    & frac_pix_sat < frac_pix_sat_thresh);  %remove saturated image cases
        
       figure(hf1)
        subplot(1,2,1); hold all; plot(Texp_us(isel)/1000, SumCount_AllPix(isel),'o-');  grid on;  
         ylabel('Sum Counts:  All Pix'); xlabel({'Texp [ms]'; '*ALL PIXELS*'}); 
        subplot(1,2,2); hold all;  plot(Texp_us(isel)/1000, SumCount_NoSatPix(isel),'o-'); grid on;  
         ylabel('Sum Counts:  NoSatPix');  xlabel({'Texp [ms]'; '* No Saturation Pixels only *'});   
        title({'Sum of Counts  v. Texp'}); 

        %***** Counts/Joule  Results *****
      figure(hf2);
        subplot(1,2,1); hold all; plot(Texp_us(isel)/1000, Count_per_J_AllPix(isel),'o-');  grid on;  
         ylabel('counts / joule:  All Pix'); xlabel({'Texp [ms]'; '*ALL PIXELS*'}); 
        subplot(1,2,2); hold all;  plot(Texp_us(isel)/1000, Count_per_J_NoSatPix(isel),'o-'); grid on;  
         ylabel('counts / joule:  NoSatPix');  xlabel({'Texp [ms]'; '* No Saturation Pixels only *'});   
        title({'Counts/Joule  v. Texp'});     
        

        %Polynomial Fit
        [params_AllPix,S] = polyfit(Texp_us(isel), Count_per_J_AllPix(isel),2);
        [params_NoSatPix,S] = polyfit(Texp_us(isel), Count_per_J_NoSatPix(isel),2);
        yfit_AllPix = polyval(params_AllPix,Texp_us(isel));
        yfit_NoSatPix = polyval(params_NoSatPix,Texp_us(isel));                        
        subplot(1,2,1);   plot(Texp_us(isel)/1000, yfit_AllPix,  '--','Color',0.7*ones(1,3));
        subplot(1,2,2);   plot(Texp_us(isel)/1000, yfit_NoSatPix,'--','Color',0.7*ones(1,3));
        %record fit values
        disp(['a[2]: ',num2str(params_AllPix(1)),'   ',num2str(params_NoSatPix(1))]);
        disp(['b[1]: ',num2str(params_AllPix(2)),'   ',num2str(params_NoSatPix(2))]);
        disp(['c[0]: ',num2str(params_AllPix(3)),'   ',num2str(params_NoSatPix(3))])
             
        fitparams.constGain.AllPix(kG,:)    = params_AllPix;
        fitparams.constGain.NoSatPix(kG,:)  = params_NoSatPix;
        fitparams.constGain.uGainVals       = uGainVals;   
        fitparams.constGain.Texp_us         = unique(Texp_us);
        
    end
    figure(hf2); legend(legnames);    figure(hf1); legend(legnames);
    %     yaxlims = ylim(); ylim([0,yaxlims(2)]);
    
    
    
 % 1D Fits v. Gain (@each Texp) ------------------------------------------------------
 
    Texp_us_rounded = round(Texp_us/1000,0)*1000;
    uTexp_us_rounded = unique( Texp_us_rounded );
    
    hf1g = figure(); hf2g = figure;
    for kT = 1:length(uTexp_us_rounded)
        Texp_sel = uTexp_us_rounded(kT); legnamesG{kT} = ['T=',num2str(Texp_sel),'us'];
        isel = find(Texp_sel == Texp_us_rounded ...
                    & frac_pix_sat < frac_pix_sat_thresh);  %remove saturated image cases

       figure(hf1g)
        subplot(1,2,1); hold all; plot(gainval(isel), SumCount_AllPix(isel),'o-');  grid on;  
         ylabel('Sum Counts:  All Pix'); xlabel({'gain'; '(*ALL PIXELS*)'})
        subplot(1,2,2); hold all;  plot(gainval(isel), SumCount_NoSatPix(isel),'o-'); grid on;  
         ylabel('Sum Counts:  NoSatPix');   xlabel({'gain'; '(Unsaturated pixels only)'});
        title({'Sum of Counts  v. Gain'; ['Max Frac. Pixels saturated allowed = ',num2str(frac_pix_sat_thresh),')']}); 

        %***** Counts/Joule  Results *****
      figure(hf2g);
        subplot(1,2,1); hold all; plot(gainval(isel), Count_per_J_AllPix(isel),'o-');  grid on;  
         ylabel('counts / joule:  All Pix');  xlabel({'gain'; '(*ALL PIXELS*)'}) 
        subplot(1,2,2); hold all;  plot(gainval(isel), Count_per_J_NoSatPix(isel),'o-'); grid on;  
         ylabel('counts / joule:  NoSatPix');  xlabel({'gain'; '(Unsaturated pixels only)'})  
        title({'Joules/Count  v. Gain'; ['(Max Frac. Pixels saturated allowed = ',num2str(frac_pix_sat_thresh),')']}); 
 
   
        %Polynomial Fit
        [params_AllPix,S] = polyfit(gainval(isel), Count_per_J_AllPix(isel),2);
        [params_NoSatPix,S] = polyfit(gainval(isel), Count_per_J_NoSatPix(isel),2);
        yfit_AllPix = polyval(params_AllPix,gainval(isel));
        yfit_NoSatPix = polyval(params_NoSatPix,gainval(isel));                        
        subplot(1,2,1);   plot(gainval(isel), yfit_AllPix,  '--','Color',0.7*ones(1,3));
        subplot(1,2,2);   plot(gainval(isel), yfit_NoSatPix,'--','Color',0.7*ones(1,3));
        %record fit values
        disp(['           ALLL PIXELS                     Un-Saturated ONLY      ']);
        disp(['a[2]: ',num2str(params_AllPix(1)),'   ',num2str(params_NoSatPix(1))]);
        disp(['b[1]: ',num2str(params_AllPix(2)),'   ',num2str(params_NoSatPix(2))]);
        disp(['c[0]: ',num2str(params_AllPix(3)),'   ',num2str(params_NoSatPix(3))])
        fitparams.constT.AllPix(kT,:) = params_AllPix;
        fitparams.constT.NoSatPix(kT,:) = params_NoSatPix;
        fitparams.constT.uTexp_us_rounded = uTexp_us_rounded;
        
        
    end
     kTT=1;
     for kT = 1: length(uTexp_us_rounded)
         legnamesG_wfits{kTT} = legnamesG{kT};
         legnamesG_wfits{kTT+1} = legnamesG{kT};  
         kTT=kTT+2;
     end 
     figure(hf2g); legend(legnamesG_wfits);    
     figure(hf1g); legend(legnamesG_wfits);
   
      
%     %***** Sum Counts  Results *****
%     figure;     plot(gainval(isel), SumCount_AllPix(isel),'o'); grid on
%     hold all;   plot(gainval(isel), SumCount_NoSatPix(isel),'o');
%     ylabel('total counts'); xlabel('gain');  legend('AllPix','NoSatPix');
%     title({'Total Counts v. Gain';['Texp[us] = ',num2str(Texp_us(1))]});
%  
%     %***** Counts / Joule  Results *****
%    
%     %Polynomial Fits
%     [params_AllPix,S] = polyfit(gainval, Count_per_J_AllPix,2);
%     [params_NoSatPix,S] = polyfit(gainval, Count_per_J_NoSatPix,2);
%     yfit_AllPix = polyval(params_AllPix,gainval);
%     yfit_NoSatPix = polyval(params_NoSatPix,gainval);
%     
%     %Counts/Joule v. Gain Plots
%     figure;     
%     plot(gainval, Count_per_J_AllPix,'o'); grid on
%     hold all;    plot(gainval, yfit_AllPix,'--');
%     hold all;   plot(gainval, Count_per_J_NoSatPix,'o');
%     hold all;     plot(gainval, yfit_NoSatPix,'--');
%     ylabel('Counts / Joule'); xlabel('gainval');  legend('AllPix','fit','NoSatPix','fit');
%     title({'Counts/Joule v. Gain';['Texp[us] = ',num2str(Texp_us(1))];['Power[W] = ',num2str(Power_W)]});
%    
%     disp(['a[2]: ',num2str(params_AllPix(1)),'   ',num2str(params_NoSatPix(1))]);
%     disp(['b[1]: ',num2str(params_AllPix(2)),'   ',num2str(params_NoSatPix(2))]);
%     disp(['c[0]: ',num2str(params_AllPix(3)),'   ',num2str(params_NoSatPix(3))]);

    %SAVE: FIT parameters
    s = split(path,'/');
    params_matrix = [fitparams.constT.NoSatPix; -111*ones(3,3);...
                     fitparams.constGain.NoSatPix; -111*ones(3,3);...
                     repmat(fitparams.twoD.params',1,3)];
    outputnameCSV= [ s{end-1},'__calibration_parameters_(NoSatPixels)_(Pin=',num2str(Power_W*1e6),'uW).csv'];
    csvwrite([dirs.input,'/',outputnameCSV],params_matrix);
    
    
    %SAVE entire workspace
    s = split(path,'/');
    outputname = [ s{end-1},'__calibration_results_(Pin=',num2str(Power_W*1e6),'uW).mat' ];
    save([dirs.input,'/',outputname]);
    
    

    % CAL FILE BUILD:  File Elements for NPZ file -- FLIR_Sensitivity.npz
    button = questdlg('generate CAL File ?')
    if strcmp(button,'Yes')
    	CAL.gain = [0:18]';
     	yfit_T1 = polyval(fitparams.constT.NoSatPix(1,:),CAL.gain); 
        yfit_T2 = polyval(fitparams.constT.NoSatPix(2,:),CAL.gain);
        %         yfit_T3 = polyval(fitparams.constT.NoSatPix(3,:),CAL.gain);
        CAL.counts_per_joule = mean([yfit_T1,yfit_T2],2);
        
      	
        %date of data collection:  from path
        expression = ['(?<month>\d+)/(?<day>\d+)/(?<year>\d+)|'...
              '(?<day>\d+)-(?<month>\d+)-(?<year>\d+)|'...
              '(?<day>\d+)_(?<month>\d+)_(?<year>\d+)'];
         token  =  regexp(path,expression)
         datestring = path(token(1):token(1)+9);
        
        %metadata defaults:
        CAL.serial            = '19415614';
        CAL.calibration_date  = datestring;
        CAL.operator          = 'Brent Fisher'
        CAL.photodiode        = 'Thorlabs S130C'
        CAL.powermeter        = 'Thorlabs PM100D ("probestation")'
        CAL.wavelength        = '940 nm'
        CAL.powermeterrange   = [num2str(730),' mW'];   
        defanswers  = {CAL.serial;CAL.calibration_date;CAL.operator;CAL.photodiode;CAL.powermeter;CAL.wavelength;CAL.powermeterrange}; 
        answers  = inputdlg({'serial';'calibdration date';'operator';'photodiode';'powermeter';'waverlength';'powermeter range'},...
                        'METADATA',1,defanswers);
        figure; plot(CAL.gain,[yfit_T1,yfit_T2]); hold all; plot(CAL.gain,CAL.counts_per_joule,'Color','k','Linewidth',3);
        xlabel('gain'); ylabel('Counts per Joule');  title('Sensitivity curve to be used in Calibration');
        legend('T=1','T=2','avg')

        s = split(path,'/');  
        outputname = ['FLIR_sensitivity_Calibration_(Pin=',num2str(Power_W*1e6),'uW)_',datestring,'.mat'];
        save([path,outputname],'CAL');
        
        disp(['Wrote file: ',outputname]);
        disp(['to:         ',path]);
        disp(CAL.counts_per_joule);
        
          
    end
    
end





function zDiff_val = zDiff_function(vars,xF,yF,zVals, fun2optimize)
    z0 = vars(1);
	ax = vars(2);
	ay = vars(3);
	bx = vars(4);
	by = vars(5);
    bxy = vars(6);
    % 	zDiff_val = sum(sum(    ( zFitpoly2_function(vars,xF,yF) - zVals ).^2    ));    
    % 	zDiff_val = sum(sum(    ( zFitpoly1_function(vars,xF,yF) - zVals ).^2    ));
%     isel = find(~isnan(xF(:)) & ~isnan(yF(:)) & ~isnan(zVals(:)) );
%     xF = xF(isel);
%     yF = yF(isel);
%     zVals = zVals(isel);
    zDiff_val = sum(sum(    ( fun2optimize(vars,xF,yF) - zVals ).^2    ));
end

function zFit = zFitpoly1_function(vars,xF,yF)
    z0 = vars(1);
	ax = vars(2);
	ay = vars(3);
    zFit = (z0 + ax*xF + ay*yF);
end

function zFit = zFitpoly2_function(vars,xF,yF)
    z0 = vars(1);
	ax = vars(2);
	ay = vars(3);
	bx = vars(4);
	by = vars(5);
    bxy = vars(6);    
    %     zFit = (z0 + ax*xF + bx*xF.^2 + ay*yF + by*yF.^2);
    zFit = (z0 + ax*xF + bx*xF.^2 + ay*yF + by*yF.^2 + bxy*xF.*yF);
end




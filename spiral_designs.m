% This plots a bunch of spiral patterns and calculates their lengths
% All Spirals start at 12:00 and rotate CLOCKWISE
% Definition of each spiral given by:
%  (a) [um] position of 12:00 on spiral (top)
%  (b) [um] Radius of Curvature of the outermost loop 
%  (c) [um] pitch 
%  (d) Number of turns  


function main()
    %                   %TOP(x,y),      ROC, pitch  , Nturns
    spirals_list{1,:}  = {[1950,  0  ] ,   300,  20,    1000 };  
    spirals_list{2,:}  = {[1400, -50 ] ,   250,  20,    1000 };  
    spirals_list{3,:}  = {[925 , -100] ,   225,  20,    1000 };  
    spirals_list{4,:}  = {[500 , -150] ,   200,  20,    1000 };  
    spirals_list{5,:}  = {[100 , -200] ,   150,  20,    1000 };  
    
    y2 = -600;
    spirals_list{6,:}  = {[1950, y2-0  ] ,   300,  12,   1000 };  
    spirals_list{7,:}  = {[1400, y2-50 ] ,   250,  12,   1000 };  
    spirals_list{8,:}  = {[925 , y2-100] ,   225,  12,   1000 };  
    spirals_list{9,:}  = {[500 , y2-150] ,   200,  12,   1000 };  
    spirals_list{10,:} = {[100 , y2-200] ,   150,  12,   1000 };  
    
    y2 = -1200;
    spirals_list{11,:}  = {[1950, y2-0  ] ,   200,  12,   1000 };  
    spirals_list{12,:}  = {[1400, y2-50 ] ,   200,  12,   1000 };  
    spirals_list{13,:}  = {[925 , y2-100] ,   200,  12,   1000 };  
    spirals_list{14,:}  = {[500 , y2-150] ,   200,  12,   1000 };  
    spirals_list{15,:}  = {[100 , y2-200] ,   200,  12,   1000 };  
 
    y2 = -1800;
    spirals_list{16,:}  = {[1950, y2-0  ] ,   250,  10,   5 };  
    spirals_list{17,:}  = {[1400, y2-50 ] ,   250,  10,   5 };  
    spirals_list{18,:}  = {[925 , y2-100] ,   250,  10,   5 };  
    spirals_list{19,:}  = {[500 , y2-150] ,   250,  10,   5 };  
    spirals_list{20,:}  = {[100 , y2-200] ,   250,  10,   5 };      
    
    
    %Loop over all defined spirals and generate a plot 
    N = size(spirals_list,1);  % how many spirals
    figure;    
    for ks = 1:N
        pos_top = spirals_list{ks}{1};
        ROC0    = spirals_list{ks}{2};   
        pitch   = spirals_list{ks}{3};   
        nturns_max = ROC0 / pitch;
        nturns 	= min([ spirals_list{ks}{4};, nturns_max]) ;
        xy_tmp = get_spiral_from_TOP_edge(pos_top,nturns, pitch, ROC0);        
        plot(xy_tmp(:,1),xy_tmp(:,2)); hold on;     
        Lspirals(ks,1) = nturns*sqrt( (2*pi*ROC0)^2 + pitch^2);       %length of spiral
        legnames{ks} = num2str(Lspirals(ks,1));
    end
    axis equal;
    legend(legnames);
    disp(Lspirals);
end

function  xy_tmp = get_spiral_from_TOP_edge(pos_top,nturns,pitch,ROC_0)
    pos_ctr(1,1) = pos_top(1,1);                 %x_ctr = x_top
    pos_ctr(1,2) = pos_top(1,2) - ROC_0;        %y_ctr = y_top - ROC_0
    phi0 = pi/2;                                %phi0 = 90deg: start at TOP
    phi = linspace(0, -nturns*2*pi, 10000) ;    % (phi<0: ClockwiseSpiral) ...  10000 = resolution
    ROC_end = max( [0, ROC_0 - nturns*pitch] );
    r = linspace(ROC_0, ROC_end, numel(phi));
    
    
    x = pos_ctr(1,1) + r .* cos(phi + phi0) ;
    y = pos_ctr(1,2) + r .* sin(phi + phi0) ;
%     plot(x,y,'b-',pos(:,1),pos(:,2),'ro-') ; % nturns crossings, including end point
    xy_tmp = [x',y'];
end


function  xy_tmp = get_spiral_from_center_to_endpoint(pos,nturns)
    dp = diff(pos,1,1) ;
    R = hypot(dp(1), dp(2)) ;
    phi0 = atan2(dp(2), dp(1)) ;
    phi = linspace(0, nturns*2*pi, 10000) ; % 10000 = resolution
    r = linspace(0, R, numel(phi)) ;
    x = pos(1,1) + r .* cos(phi + phi0) ;
    y = pos(1,2) + r .* sin(phi + phi0) ;
%     plot(x,y,'b-',pos(:,1),pos(:,2),'ro-') ; % nturns crossings, including end point
    xy_tmp = [x',y'];
end


 




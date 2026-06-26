% generate scttering distribution for waveguide sidewall
%

function generate_waveguide_ies()
    % --- PHYSICAL PARAMETERS ---
    n_eff = 2.25;          % TM04 mode effective index
    n_clad = 1.45;          % Air cladding
    exponent_n = 6;        % Sharpness of the scattering cone (1=Lambertian)
    
    % Calculate the forward-bias peak angle (Cherenkov/Phase condition)
    % Measured relative to the waveguide axis (Z).
    theta_peak_rad = acos(n_clad / n_eff); 
    theta_peak_deg = rad2deg(theta_peak_rad); 
    
    % The target tilt angle for the IES grid
    target_tilt = theta_peak_deg; 

    % --- DYNAMIC FILENAME GENERATION ---
    filename = sprintf('waveguide_scat_neff%.2f_n%d.ies', n_eff, exponent_n);

    % --- ANGULAR GRID DEFINITION ---
    vert_angles = 0:1:180;      % 1-degree steps for a smoother plot/file
    horiz_angles = 0:10:360;
    
    num_vert = length(vert_angles);
    num_horiz = length(horiz_angles);
    
    % --- CALCULATE INTENSITY MATRIX (Global Z = Forward, Global X = Normal) ---
    candela_matrix = zeros(num_vert, num_horiz);
    
    % Target angle measured from the normal (Local Z) toward the forward axis (Local X)
    theta_normal_target = 90 - theta_peak_deg; 
    
    % Direction vector of the ideal peak scattering lobe in IES local space
    tx = sind(theta_normal_target);
    ty = 0;
    tz = cosd(theta_normal_target);
    
    for h = 1:num_horiz
        phi = horiz_angles(h);
        for v = 1:num_vert
            theta = vert_angles(v);
            
            % Convert standard IES spherical grid to a local 3D unit direction vector
            kx = sind(theta) * cosd(phi); % Local X (Maps to Global +Z Forward)
            ky = sind(theta) * sind(phi); % Local Y (Maps to Global -Y)
            kz = cosd(theta);             % Local Z (Maps to Global +X Normal)
            
            % 1. Hemispheric Mask: Light can only emerge out into the air (kz > 0)
            if kz > 0
                % 2. Calculate cosmic distance (dot product) from our target peak lobe
                cos_dist = (kx*tx + ky*ty + kz*tz);
                
                if cos_dist > 0
                    % Apply the scattering cone profile
                    candela_matrix(v, h) = cos_dist^exponent_n;
                else
                    candela_matrix(v, h) = 0;
                end
            else
                candela_matrix(v, h) = 0; % Spliced out (trapped inside waveguide)
            end
        end
    end
    
    % Normalize to a peak of 1000 candelas
    candela_matrix = (candela_matrix / max(candela_matrix(:))) * 1000;

    % --- WRITE IESNA FILE FORMAT ---
    fid = fopen(filename, 'w');
    fprintf(fid, 'IESNA:LM-63-2002\n');
    fprintf(fid, '[TEST] Custom Waveguide Sidewall Scattering\n');
    fprintf(fid, '[PARAM] n_eff=%.2f, exponent=%d, peak_angle=%.1f deg\n', n_eff, exponent_n, target_tilt);
    fprintf(fid, 'TILT=NONE\n');
    fprintf(fid, '1 -1 1 %d %d 1 1 0 0 0\n', num_vert, num_horiz); 
    fprintf(fid, '1.0 1.0 0.0\n'); 
    
    fprintf(fid, '%f ', vert_angles); fprintf(fid, '\n');
    fprintf(fid, '%f ', horiz_angles); fprintf(fid, '\n');
    
    for h = 1:num_horiz
        fprintf(fid, '%f ', candela_matrix(:, h));
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('Successfully generated file: %s\n', filename);

% % --- GENERATE POLAR PLOT VISUALIZATION ---
%     % Grab a slice of the distribution along a single horizontal slice (phi = 0)
%     plot_angles_rad = deg2rad(vert_angles);
%     intensity_slice = candela_matrix(:, 1); 
%     
%     figure('Color', 'w', 'Name', 'IESNA Far-Field Distribution');
%     polarplot(plot_angles_rad, intensity_slice, 'r-', 'LineWidth', 2);
%     
%     % Format the polar plot so 0 degrees is the forward waveguide axis
%     ax = gca;
%     ax.ThetaZeroLocation = 'top'; % Set 0 degrees at the top vertical axis
%     ax.ThetaDir = 'clockwise';
%     ax.ThetaLim = [0 180];        % Show forward hemisphere
%     
%     % Fixed newline interpretation using a cell array
%     title_str = {sprintf('Scattering Distribution Profile'), ...
%                  sprintf('(n_{eff} = %.2f, Peak Lobe = %.1f^\\circ)', n_eff, target_tilt)};
%     title(title_str, 'FontSize', 12);
%     
%     % Correct way to guide the user on the radial axis units
%     text(deg2rad(15), max(intensity_slice)*0.6, 'Intensity (Candela)', ...
%          'Color', [0.4 0.4 0.4], 'FontWeight', 'bold');
% --- GENERATE 2D CONTOUR PLOT VISUALIZATION ---
    figure('Color', 'w', 'Name', '2D Far-Field Mapping', 'Position', [150 150 700 550]);
    
    % Create a 2D grid of the horizontal and vertical angles for plotting
    [Horiz_Grid, Vert_Grid] = meshgrid(horiz_angles, vert_angles);
    
    % Plot using filled contours
    contourf(Horiz_Grid, Vert_Grid, candela_matrix, 20, 'LineColor', 'none');
    colormap('hot'); % 'hot' or 'jet' works beautifully for intensity maps
    colorbar;
    
    % Format the axes to match the physical meaning
    set(gca, 'YDir', 'reverse'); % Traditional IES viewing: 0 at top, 180 at bottom
    xlabel('Horizontal Angle (\phi) [Degrees]', 'FontWeight', 'bold');
    ylabel('Vertical Angle (\theta) [Degrees]', 'FontWeight', 'bold');
    
    % Add custom markings to make the coordinate system obvious
    hold on;
    % Draw a line at the peak vertical angle
    yline(target_tilt, 'w--', sprintf('Peak Lobe (\\theta = %.1f^\\circ)', target_tilt), ...
          'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.8 0.8 0.8]);
      
    % Highlight the boundaries of the allowed outward hemisphere (\phi <= 90 or \phi >= 270)
    xline(90, 'w:', 'Outward Boundary', 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7]);
    xline(270, 'w:', 'Outward Boundary', 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7]);
    
    title_str_2d = {sprintf('2D Source Intensity Mapping (Candela Scale)'), ...
                    sprintf('Forward Vector = [\\theta=%.1f^\\circ, \\phi=0^\\circ/360^\\circ]', target_tilt)};
    title(title_str_2d, 'FontSize', 12);
    
    axis tight;


end

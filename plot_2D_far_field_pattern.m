
function  plot_2D_far_field_pattern()
    %% 1. Path and File Selection
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_03_02_LabData\nichia_NA';

    % Check if directory exists before changing
    if exist(dirs.d, 'dir')
        cd(dirs.d);
    else
        fprintf('Path not found. Please select folder manually.\n');
        dirs.d = uigetdir();
        cd(dirs.d);
    end

    [files, path] = uigetfile('*.txt', 'Select X-axis then Y-axis file', 'Multiselect', 'on');

    % Handle cases where user cancels or forgets to select two files
    if isequal(files,0) || length(files) < 2
        error('Please select at least two files (Theta_X and Theta_Y).');
    end

    %% 2. Load Data (Using readtable for R2018b compatibility)
    % readtable is more robust than csvread for files with headers
    data_x_table = readtable(fullfile(path, files{1}));
    data_y_table = readtable(fullfile(path, files{2}));

    % Convert table to array (Assuming Col 1 is Angle, Col 2 is Intensity)
    data_x = table2array(data_x_table);
    data_y = table2array(data_y_table);

    theta_x_deg = data_x(:,1);
    int_x = data_x(:,2);

    theta_y_deg = data_y(:,1);
    int_y = data_y(:,2);

    %% 3. Create 2D Intensity Map & NA Conversion
    % Using outer product to create the elliptical 2D distribution
    Intensity_Map = int_y * int_x'; 

    % Convert to Numerical Aperture: NA = sin(theta)
    na_x = sin(deg2rad(theta_x_deg));
    na_y = sin(deg2rad(theta_y_deg));
    [NA_X, NA_Y] = meshgrid(na_x, na_y);

    % Calculate radial NA distance for each pixel
    NA_R_Map = sqrt(NA_X.^2 + NA_Y.^2);

    %% 4. Integrate Encircled Energy
    na_r_axis = linspace(0, max([max(na_x), max(na_y)]), 100);
    encircled_fraction = zeros(size(na_r_axis));
    total_intensity = sum(Intensity_Map(:));

    for i = 1:length(na_r_axis)
        mask = NA_R_Map <= na_r_axis(i);
        encircled_fraction(i) = sum(Intensity_Map(mask)) / total_intensity;
    end

    %% 5. Visualization
    figure('Color', 'w', 'Position', [100, 100, 1100, 450]);

    % Plot 2D Far Field
    subplot(1,2,1);
    imagesc(na_x, na_y, Intensity_Map);
    axis image; colorbar;
    colormap('jet');
    xlabel('NA_x'); ylabel('NA_y');
    title('2D Far-Field Intensity (NA Space)');

    % Plot Encircled Energy Curve
    subplot(1,2,2);
    plot(na_r_axis, encircled_fraction, 'b-', 'LineWidth', 2);
    grid on; hold on;
    xlabel('NA Radius (Size of Aperture)');
    ylabel('Fraction of Enclosed Power');
    title('Encircled Energy vs. NA');


end
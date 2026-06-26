
function compare_defect_lists()
    %% Setup and Parameters
    root_dir = 'C:\Users\brent\Downloads\defect_compare';
    x_tolerance = 50.0;      
    mag_tolerance = 100.0;    
    
    % 1. Select the first file
    cd(root_dir);
    [f1, p1] = uigetfile('*.csv', 'Select CSV for Sample 1');
    if isequal(f1,0), error('User cancelled file selection for Sample 1'); end
    path1 = fullfile(p1, f1);

    % 2. Select the second file
    cd(root_dir);
    [f2, p2] = uigetfile('*.csv', 'Select CSV for Sample 2');
    if isequal(f2,0), error('User cancelled file selection for Sample 2'); end
    path2 = fullfile(p2, f2);

    % 3. Load the data (R2018b compatible syntax)
    % In 2018b, we use 'PreserveVariableNames' instead of 'VariableNamingRule'
    opts1 = detectImportOptions(path1);
    T1 = readtable(path1, opts1);

    opts2 = detectImportOptions(path2);
    T2 = readtable(path2, opts2);

%% Comparison Logic
    all_ids = unique([T1.("Waveguide"); T2.("Waveguide")]);
    results = table('Size', [length(all_ids), 2], ...
        'VariableTypes', {'string', 'double'}, ...
        'VariableNames', {'WaveguideID', 'MatchFraction'});

    for i = 1:length(all_ids)
        current_id = all_ids{i};
        rows1 = T1(strcmp(T1.("Waveguide"), current_id), :);
        rows2 = T2(strcmp(T2.("Waveguide"), current_id), :);

        n1 = height(rows1); n2 = height(rows2);

        if n1 == 0 || n2 == 0
            results.WaveguideID(i) = current_id;
            results.MatchFraction(i) = 0;
            continue;
        end

        match_count = 0;
        sum_matched_mags = 0;
        matched_in_2 = false(n2, 1);

        for r1 = 1:n1
            for r2 = 1:n2
                if matched_in_2(r2), continue; end

                if (abs(rows1.x_um(r1) - rows2.x_um(r2)) <= x_tolerance) && ...
                   (abs(rows1.mag(r1) - rows2.mag(r2)) <= mag_tolerance)

                    % Standard Match
                    match_count = match_count + 1;

                    % Weighted Match: Use mean magnitude of the pair
                    mean_mag = (rows1.mag(r1) + rows2.mag(r2)) / 2;
                    sum_matched_mags = sum_matched_mags + mean_mag;

                    matched_in_2(r2) = true;
                    break;
                end
            end
        end        
        
        % Calculation Logic:
        % Standard Fraction: Matches / Max potential number of defects
        results.MatchFraction(i) = match_count / max(n1, n2);

        % Weighted Fraction: Total matched energy / Max potential energy
        total_mag_1 = sum(rows1.mag);
        total_mag_2 = sum(rows2.mag);
        max_possible_mag = max(total_mag_1, total_mag_2);

        results.WeightedMatch(i) = sum_matched_mags / max_possible_mag;
        results.WaveguideID(i) = current_id;
        
        
    end

    %% Visualization Section
    figure('Color', 'w', 'Position', [100, 100, 1000, 500]);

    % --- Subplot 1: Bar Chart of Match Fractions ---
    subplot(1, 2, 1);
    hBar = bar(categorical(results.WaveguideID), results.MatchFraction);
    hBar.FaceColor = 'flat';
    hBar.CData = repmat([0.2 0.6 0.8], height(results), 1); % Nice blue color

    ylim([0 1.1]); % Extra room for labels
    grid on;
    xlabel('Waveguide ID');
    ylabel('Match Fraction (0.0 to 1.0)');
    title('Defect Match Consistency per Waveguide');

    % --- Subplot 2: Comparison Scatter (Visualizing X and Mag) ---
    subplot(1, 2, 2);
    hold on;
    scatter(T1.x_um, T1.mag, 60, 'blue', 'filled', 'MarkerFaceAlpha', 0.5);
    scatter(T2.x_um, T2.mag, 40, 'red', 'DisplayName', 'Sample 2');
    set(gca, 'YScale', 'log')
    hold off;

    grid on;
    xlabel('X Location (\mum)');
    ylabel('Magnitude');
    title('Overlay of All Defect Locations');
    legend('Sample 1', 'Sample 2', 'Location', 'best');

    % Adjust layout
    sgtitle({['Chip Comparison (Tol_x: ' num2str(x_tolerance) 'um,  Tol_m: ' num2str(mag_tolerance) ')'];...
        ['file1: ',replace(f1,'_',' ')];['file2: ',replace(f2,'_',' ')]});
    
        %% Visualization Section 2
    figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
    hBar = bar(categorical(results.WaveguideID), results.WeightedMatch);
    hBar.FaceColor = 'flat';
    hBar.CData = repmat([0.2 0.6 0.8], height(results), 1); % Nice blue color
    ylim([0 1.1]); % Extra room for labels
    grid on;
    xlabel('Waveguide ID');
    ylabel('Match Fraction (magnitude-weighted; 0.0 to 1.0)');
    title({'Defect Match Consistency per Waveguide'; ...
        ['Chip Comparison (Tol_x: ' num2str(x_tolerance) 'um,  Tol_m: ' num2str(mag_tolerance) ')'];...
        ['file1: ',replace(f1,'_',' ')];['file2: ',replace(f2,'_',' ')]});
    
 
end
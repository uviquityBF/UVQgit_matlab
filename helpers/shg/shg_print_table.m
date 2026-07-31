function shg_print_table(cases, blue_losses, matrix_W, num_cases, num_losses)
% Prints a Configuration x Loss table of SHG power values (W) as microwatts.
    fprintf('%-36s', 'Configuration');
    for c = 1:num_losses
        fprintf(' | a0=%2ddB/cm', blue_losses(c));
    end
    fprintf('\n%s\n', repmat('-', 38 + num_losses * 13, 1));
    for k = 1:num_cases
        fprintf('%-36s', sprintf('%d: %s', k, cases(k).name));
        for c = 1:num_losses
            fprintf(' | %10.3f uW', matrix_W(k, c) * 1e6);
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-', 38 + num_losses * 13, 1));
end

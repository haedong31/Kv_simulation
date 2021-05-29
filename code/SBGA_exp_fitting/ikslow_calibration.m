function [par] = ikslow_calibration(amp, tau, tol, N0, N1, N2)
    % calibration arguments
    global num_var;
    global low;
    global high;
    global max_iter;
    global window_size;

    global hold_volt;
    global hold_t;
    global volts;
    global pulse_t;
    global Ek;

    num_var = 6;
    low = [-90.0, 2.0, -20.0, 2.0, 200.0, 0.01]; 
    high = [90.0, 14.0, 80.0, 24.0, 5000.0, 0.31];
    max_iter = 100;
    window_size = N1; % N1 as pooling window size

    hold_volt = -70;
    hold_t = 0.125*1000;
    volts = 0:10:50;
    pulse_t = 4.5*1000;
    Ek = -91.1;

    % run calibration
    par = kslow_sbga(amp, tau, tol, N0, N1, N2);
end

function [best_chroms] = kslow_sbga(amp, tau, tol, N0, N1, N2)
    global num_var;
    global max_iter;
    global window_size;

    % initial run
    iter = 1;
    chroms = init_gen(N0);
    sig_list = zeros(window_size, num_var);
    eval_list = zeros(max_iter, 1);

    % evaluation
    [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
    total_diff = amp_diff + tau_diff;
    [min_diff, min_diff_idx] = min(total_diff);
    eval_list(iter) = min_diff;
    fprintf('Initial best fit: %f|Amp: %f|Tau: %f \n', min_diff, amp_diff(min_diff_idx), tau_diff(min_diff_idx));

    % evolution
    [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);

    % repeat; tolerance or max_iter
    while true
        if iter > max_iter
            iter = 1;
            chroms = init_gen(N0);
            sig_list = zeros(window_size, num_var);
            eval_list = zeros(max_iter, 1);
        
            % evaluation
            [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
            total_diff = amp_diff + tau_diff;
            [min_diff, min_diff_idx] = min(total_diff);
            eval_list(iter) = min_diff;
            fprintf('Initial best fit: %f|Amp: %f|Tau: %f \n', min_diff, amp_diff(min_diff_idx), tau_diff(min_diff_idx));
        
            % evolution
            [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);
        end

        iter = iter + 1;
        fprintf('Generation %i \n', iter)

        % evaluation
        [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
        total_diff = amp_diff + tau_diff;
        [min_diff, min_diff_idx] = min(total_diff);
        eval_list(iter) = min_diff;        

        % check stopping tolerance
        min_amp_diff = amp_diff(min_diff_idx);
        min_tau_diff = tau_diff(min_diff_idx);
        if ((min_amp_diff <= tol(1)) && (min_tau_diff <= tol(2)))
            fprintf('Termination: %f|Amp: %f|Tau: %f \n', min_diff, min_amp_diff, min_tau_diff);
            best_chroms = chroms(min_diff_idx, :);
            break
        end

        % check solution has been improved
        if (min_diff < eval_list(iter-1))
            fprintf('Updated best fit: %f|Amp: %f|Tau: %f \n', min_diff, min_amp_diff, min_tau_diff);
        end

        % evolution
        [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);
    end
end

function [chroms] = init_gen(N0)
    global num_var;
    global low;
    global high;

    chroms = zeros(N0, num_var);
    for j = 1:num_var
        unif_dist = makedist('Uniform', 'lower',low(j), 'upper',high(j));
        chroms(:, j) = random(unif_dist, N0, 1);
    end
end

function [new_chroms, new_sig_list] = next_gen(chroms, evals, sig_list, N0, N1, N2)
    global num_var;
    global window_size

    % superior chromosomes 
    [~, sup_chrom_idx] = mink(evals, N1);
    sup_chroms = chroms(sup_chrom_idx, :);

    % renew sigma list
    new_sig_list = sig_list; % copy
    sig = std(sup_chroms);
    new_sig_list(1:(window_size-1), :) = new_sig_list(2:window_size, :);
    new_sig_list(window_size, :) = sig;

    % pooled sigma
    pooled_sig = mean(new_sig_list);

    % mean chromosome of superiors
    mean_elite = mean(sup_chroms, 1);
    sup_chroms(end, :) = mean_elite;

    % include superiors in new generation
    new_chroms = zeros(N0, num_var);
    new_chroms(1:N1, :) = sup_chroms;
    
    % breed new chromosomes from superiors
    cnt = 1;
    for i = 1:N1
        sup_chrom = sup_chroms(i, :);
        for j = 1:N2
            new_chrom = sup_chrom + normrnd(0, pooled_sig);
            new_chroms((N1+cnt), :) = new_chrom;
            cnt = cnt + 1;
        end
    end
end

function [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0)
    global hold_volt;
    global hold_t;
    global volts;
    global pulse_t;
    global Ek;

    num_volts = length(volts);
    amp_diff = zeros(N0, num_volts);
    tau_diff = zeros(N0, num_volts);
    
    for i = 1:N0
        for j = 1:num_volts
            volt = volts(j);
            try
                [t, ~, A] = ikto(chroms(i, :), hold_volt, hold_t, volt, pulse_t, Ek);
                current_trace = A(:, 5);
                
                % check validity of trace shape
                [peak, peak_idx] = max(current_trace);
                [~, hold_idx] = min(abs(t - hold_t));

                check_pt1 = any(isnan(current_trace));
                check_pt2 = any(current_trace < 0); 
                check_pt3 = var(current_trace(1:hold_idx)) > 1; % not stable at hold_volt
                check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

                if (check_pt1 || check_pt2 || check_pt3 || check_pt4)
                    amp_diff(i, j) = 1e+4;
                    tau_diff(i, j) = 1e+4;
                else
                    current_trace_trunc = current_trace(peak_idx:end);
                    t_trunc = t(peak_idx:end);
                    t_trunc = t_trunc - t_trunc(1);

                    [~, tau_idx] = min(abs(peak*exp(-1) - current_trace_trunc));
                    tau_hat = t_trunc(tau_idx);

                    % calculate discrepancies
                    amp_diff(i, j) = abs(peak - amp(j));
                    tau_diff(i, j) = abs(tau_hat - tau(j));
                end
            catch
                % large discrepancy for invalid trace shape
                amp_diff(i, j) = 1e+4;
                tau_diff(i, j) = 1e+4;
            end
        end
    end
    amp_diff = mean(amp_diff, 2);
    tau_diff = mean(tau_diff, 2);
end

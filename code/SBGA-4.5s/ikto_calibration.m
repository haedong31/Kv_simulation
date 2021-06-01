function [par, amp_diff, tau_diff] = ikto_calibration(amp, tau, input_t, input_volt, N1, N2)
    % calibration arguments
    global num_var;
    global low;
    global high;
    global max_iter;
    global log_interval;
    global window_size;

    global t;
    global hold_volt;
    global hold_idx;
    global volt;
    global Ek;

    num_var = 5;
    low = [0, 0, 0, 1, 0];
    high = [60, 50, 60, 15, 0.5];
    max_iter = 300;
    log_interval = 50;
    window_size = N1; % N1 as pooling window size

    t = input_t;
    hold_volt = -70;
    volt = input_volt;
    Ek = -91.1;

    % estimate index of holding time
    ideal_hold_time = 0.125*1000;
    [~, hold_idx] = min(abs(t - ideal_hold_time));    

    % run calibration
    N0 = N1 + N1*N2;
    [par, amp_diff, tau_diff] = kto_sbga(amp, tau, N0, N1, N2);
end

function [last_chrom, last_amp_diff, last_tau_diff] = kto_sbga(amp, tau, N0, N1, N2)
    global num_var;
    global max_iter;
    global window_size;
    global log_interval;

    % initial run
    iter = 1;
    chroms = init_gen(N0);
    sig_list = zeros(window_size, num_var);
    
    % evaluation
    [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
    total_diff = amp_diff + tau_diff;
    [min_diff, min_diff_idx] = min(total_diff);
    min_amp_diff = amp_diff(min_diff_idx);
    min_tau_diff = tau_diff(min_diff_idx);   

    % print evaluation results
    fprintf('Initial fit: %f|Amp: %f|Tau: %f \n', min_diff, min_amp_diff, min_tau_diff);

    % evolution
    [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);

    % repeat; tolerance or max_iter
    while iter <= max_iter
        % increase interation
        iter = iter + 1;

        % evaluation
        [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
        total_diff = amp_diff + tau_diff;
        
        % check point; print discrepancies
        if (mod(iter, log_interval) == 0)
            [min_diff, min_diff_idx] = min(total_diff);
            min_amp_diff = amp_diff(min_diff_idx);
            min_tau_diff = tau_diff(min_diff_idx);    

            fprintf('[%i/%i] %f|Amp: %f|Tau: %f \n', iter, max_iter, min_diff, min_amp_diff, min_tau_diff);
        end
        
        % evolution
        [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);
    end
    
    [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
    total_diff = amp_diff + tau_diff;
    [~, min_diff_idx] = min(total_diff);
            
    last_chrom = chroms(min_diff_idx, :);
    last_amp_diff = amp_diff(min_diff_idx);
    last_tau_diff = tau_diff(min_diff_idx); 
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
    global t;
    global hold_volt;
    global hold_idx;
    global volt;
    global Ek;

    hold_t = t(1:hold_idx);
    pulse_t = t((hold_idx+1):end);
    pulse_t_adj = pulse_t - pulse_t(1);
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = hold_t;
    time_space{3} = pulse_t_adj;

    amp_diff = zeros(N0, 1);
    tau_diff = zeros(N0, 1);
    
    for i = 1:N0
        try
            current_trace = ikto(chroms(i, :), hold_volt, volt, time_space, Ek);
            
            % check validity of trace shape
            [peak, peak_idx] = max(current_trace);

            check_pt1 = any(isnan(current_trace));
            check_pt2 = any(current_trace < 0); 
            check_pt3 = var(current_trace(1:hold_idx)) > 0.01; % not stable at hold_volt
            check_pt4 = peak_idx < hold_idx; % not stable at hold_volt of too flat at pulse

            if (check_pt1 || check_pt2 || check_pt3 || check_pt4)
                amp_diff(i) = 1e+4;
                tau_diff(i) = 1e+4;
            else
                current_trace_trunc = current_trace(peak_idx:end);
                t_trunc = t(peak_idx:end);
                t_trunc = t_trunc - t_trunc(1);

                [~, tau_idx] = min(abs(peak*exp(-1) - current_trace_trunc));
                tau_hat = t_trunc(tau_idx);

                % calculate discrepancies
                amp_diff(i) = abs(peak - amp);
                tau_diff(i) = abs(tau_hat - tau);
            end
        catch
            % large discrepancy for invalid trace shape
            amp_diff(i) = 1e+4;
            tau_diff(i) = 1e+4;
        end
    end
end

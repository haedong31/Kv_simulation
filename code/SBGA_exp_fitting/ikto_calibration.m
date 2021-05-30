function [par] = ikto_calibration(amp, tau, volts_input, N1, N2, verbose)
    % calibration arguments
    global num_var;
    global low;
    global high;
    global max_iter;
    global log_interval;
    global window_size;

    global hold_volt;
    global hold_len;
    global volts;
    global end_len;
    global Ek;

    num_var = 5;
    low = [-40, 0, -40, 1, 0];
    high = [60, 50, 60, 15, 0.5];
    max_iter = 500;
    log_interval = 50;
    window_size = N1; % N1 as pooling window size

    hold_volt = -70;
    hold_len = 0.125*1000;
    volts = volts_input;
    end_len = 4.5*1000;
    Ek = -91.1;

    % run calibration
    N0 = N1 + N1*N2;
    par = kto_sbga(amp, tau, N0, N1, N2, verbose);
end

function [best_chrom] = kto_sbga(amp, tau, N0, N1, N2, verbose)
    global num_var;
    global max_iter;
    global window_size;
    global log_interval;

    % initial run
    iter = 1;
    chroms = init_gen(N0);
    sig_list = zeros(window_size, num_var);

    total_diff_list = [];
    amp_diff_list = [];
    tau_diff_list  = [];
    
    % evaluation
    [amp_diff, tau_diff] = evaluation(amp, tau, chroms, N0);
    total_diff = amp_diff + tau_diff;
    [min_diff, min_diff_idx] = min(total_diff);
    min_amp_diff = amp_diff(min_diff_idx);
    min_tau_diff = tau_diff(min_diff_idx);   

    % print evaluation results
    total_diff_list = [total_diff_list; min_diff];
    amp_diff_list = [amp_diff_list; min_amp_diff];
    tau_diff_list = [tau_diff_list; min_tau_diff];
    fprintf('Initial fit: %f|Amp: %f|Tau: %f \n', min_diff, min_amp_diff, min_tau_diff);

    % draw figure and linkdata
    if verbose == true
        figure(1)
        plot(total_diff_list);
        title('Total')
        
        figure(2)
        plot(amp_diff_list);
        title('Amp')
        
        figure(3)
        plot(tau_diff_list);
        title('Tau')
        
        linkdata on
    end

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

            total_diff_list = [total_diff_list; min_diff];
            amp_diff_list = [amp_diff_list; min_amp_diff];
            tau_diff_list = [tau_diff_list; min_tau_diff];
            fprintf('[%i/%i] %f|Amp: %f|Tau: %f \n', iter, max_iter, min_diff, min_amp_diff, min_tau_diff);

            if verbose == true
                figure(1)
                plot(total_diff_list);
                title('Total')
            
                figure(2)
                plot(amp_diff_list);
                title('Amp')
            
                figure(3)
                plot(tau_diff_list);
                title('Tau')
            end
        end
        
        % evolution
        [chroms, sig_list] = next_gen(chroms, total_diff, sig_list, N0, N1, N2);
    end
    
    best_chrom = chroms(min_diff_idx, :);
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
    global hold_len;
    global volts;
    global end_len;
    global Ek;

    hold_t = 0:hold_len;
    pulse_t = (hold_len + 1):end_len;
    pulse_t_adj = pulse_t - pulse_t(1);
    t = [hold_t, pulse_t];
    
    time_space = cell(1,3);
    time_space{1} = t;
    time_space{2} = hold_t;
    time_space{3} = pulse_t_adj;

    num_volts = length(volts);
    amp_diff = zeros(N0, num_volts);
    tau_diff = zeros(N0, num_volts);
    
    for i = 1:N0
        for j = 1:num_volts
            volt = volts(j);
            try
                current_trace = ikto(chroms(i, :), hold_volt, volt, time_space, Ek);
                
                % check validity of trace shape
                [peak, peak_idx] = max(current_trace);
                [~, hold_idx] = min(abs(t - hold_len));

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

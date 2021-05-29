function [par] = ito_calibration(t, yksum, tol, N0, N1, N2)
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
    low = [0.0, 0.0, 0.0, 20.0, 2.0, 0.2];
    high = [70.0, 70.0, 70.0, 70.0, 50.0, 0.6];
    max_iter = 100;
    window_size = N1; % N1 as pooling window size

    hold_volt = -70;
    hold_t = 0.125*1000;
    volts = 0:10:50;
    pulse_t = 4.5*1000;
    Ek = -91.1;

    num_volts = length(volts);
    amp = zeros(3, num_volts);
    tau = zeros(2, num_volts);

    for i = 1:num_volts
        [amp_running, tau_running] = bi_exp_fit(t, yksum(:, i));
        amp(:, i) = amp_running;
        tau(:, i) = tau_running;
    end

    % run calibration
    best_chroms = ito_sbga(amp, tau, tol, N0, N1, N2);
    par = best_chroms(end, :);
end

function [best_chroms] = ito_sbga(amp, tau, tol, N0, N1, N2)
    global num_var;
    global low;
    global high;
    global max_iter;
    global window_size;

    sig_list = zeros(window_size, num_var);

end

function [chroms] = init_gen(N0)
    global num_var;
    global low;
    global high;

    chroms = zeros(N0, num_var);
    for j = 1:num_var
        unif_dist = makedist('Uniform', 'lower',low(j) 'upper');
        chroms(:, j) = random(unif_dist, N0, 1);
    end
end

function [new_chroms] = evolve(chroms, evals, sig_list, N0, N1, N2)
    global num_var;

    % superior chromosomes 
    [~, sup_chrom_idx] = mink(evals, N1);
    sup_chroms = chroms(sup_chrom_idx, :);

    % pooled sigma
    sig = std(sup_chroms);
    pooled_sig = mean(sig_list);

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
                [t, ~, A] = Ito(chroms(i, :), hold_volt, hold_t, volt, pulse_t, Ek);
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
                    % calculate discrepancies
                    amp_diff(i, j) = abs(peak - )
                    tau_diff(i, j) = 
                end
            catch
                % large discrepancy for invalid trace shape
                amp_diff(i, j) = 1e+4;
                tau_diff(i, j) = 1e+4;
            end
        end
    end
end

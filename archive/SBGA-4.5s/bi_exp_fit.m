function [amp, tau] = bi_exp_fit(t, yksum)
    % cut early increasing phase
    [peak, peak_idx] = max(yksum);
    t_trunc = t(peak_idx:end);
    t_trunc = t_trunc - t_trunc(1);
    yksum_trunc = yksum(peak_idx:end);
    
    % objective function
    obj = @(x) exp_rmse_obj(x, yksum_trunc, t_trunc);
    
    % initial value
    x0 = rand(5, 1);

    % initial value for amplitudes; split of peak of yksum
    x0(1:3) = x0(1:3)/sum(x0(1:3));
    x0(1:3) = peak*x0(1:3);

    % initial value for tau; scaled random number by tau of yksum
    [~, tau_idx] = min(abs(yksum_trunc - peak*exp(-1)));
    tau = t_trunc(tau_idx);
    x0(4:5) = tau*x0(4:5);

    % constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];

    % lower and upper bound constraints
    lb = ones(5,1)*eps;
    ub = [ones(3,1)*peak; ones(2,1)*t_trunc(end)];

    % run optimization
    options = optimoptions('fmincon', 'Display','off');
    [sol, ~] = fmincon(obj, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

    amp = sol(1:3);
    tau = sol(4:5);
end

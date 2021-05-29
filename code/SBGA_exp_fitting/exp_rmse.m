function [rmse] = exp_rmse(x, y, t)
    amp_to = x(1);
    amp_kslow = x(2);
    amp_ss = x(3);
    tau_to = x(4);
    tau_kslow = x(5);

    yhat = amp_to.*exp(-t./tau_to) + amp_kslow.*exp(-t./tau_kslow) + amp_ss;
    rmse = sqrt(mean((y - yhat).^2));
end

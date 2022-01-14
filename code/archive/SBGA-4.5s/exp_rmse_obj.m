function [rmse] = exp_rmse_obj(x, y, t)
    amp1 = x(1);
    amp2 = x(2);
    amp_ss = x(3);
    tau1 = x(4);
    tau2 = x(5);

    yhat = amp1.*exp(-t./tau1) + amp2.*exp(-t./tau2) + amp_ss;
    rmse = sqrt(mean((y - yhat).^2));
end

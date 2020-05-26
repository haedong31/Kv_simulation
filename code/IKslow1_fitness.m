function e = IKslow1_fitness(X, Y, t, V)
    IKslow1_trc = IKslow1(t, V, X);

    IKslow1_hat = max(IKslow1_trc);
    amp_slow1 = nansum((Y.AMP - IKslow1_hat).^2);

    [~, tau1_idx] = min(abs(IKslow1_hat*exp(-1) - IKslow1_trc));
    tau_slow1 = nansum((Y.TAU - t(tau1_idx)).^2);
    
    e = amp_slow1 + tau_slow1;
end

function z = Ktrace_fitness(X, IKslow1, tau1)    
    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    Ek = -80.3;

    [t, ~, A, ~] = IKslow(X, holding_p, holding_t, P1, P1_t, Ek);
    IKslow1_trc = A(:,5);
    
    [peak, ~] = max(IKslow1_trc);
    [~, tau_idx] = min(abs(peak*exp(-1) - IKslow1_trc));
    tau = t(tau_idx);
    
    amp_delta = (peak - IKslow1).^2;
    tau_delta = (tau - tau1).^2;
    e = amp_delta + tau_delta
    fprintf('Amp: %f Tau: %f SSE: %f \n', amp_delta, tau_delta, e)
end   

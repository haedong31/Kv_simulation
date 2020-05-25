function z = Ktrace_fitness(X, ds_Ktrace, Iss, Ito, IKslow1, IKslow2, tau_to, tau1, tau2)
    % protocol
    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    P2 = -70; % mV
    P2_t = P1_t; % ms

    [t, ~, A, ~] = Kv(X, holding_p, holding_t, P1, P1_t, P2, P2_t);

    % IKsum without Iss
    Ito_trc = A(:,13);
    IKslow1_trc = A(:,14);
    IKslow2_trc = A(:,15);
    IKsum = Ito_trc + IKslow1_trc + IKslow2_trc;

    % add the constant term (Iss) after the peak
    [~, peak_idx] = max(IKsum);
    IKsum(peak_idx:end) = IKsum(peak_idx:end) + Iss;

    % align the two traces and downsampling IKsum
%     IKsum = alignsignals(IKsum, ds_Ktrace.I);
%     IKsum_tbl = table(t, IKsum);
%     sample_rate = floor(length(t)/length(ds_Ktrace.time));
% 
%     if sample_rate > 0
%         ds_IKsum = downsample(IKsum_tbl, sample_rate);
%     else
%         ds_IKsum = IKsum_tbl;
%     end

    trace_sim = dtw(IKsum, ds_Ktrace.I);

    Ito_hat = max(Ito_trc);
    IKslow1_hat = max(IKslow1_trc);
    IKslow2_hat = max(IKslow2_trc);

    amp_to = abs(Ito - max(Ito_trc));
    amp_slow1 = abs(IKslow1 - max(IKslow1_trc));
    amp_slow2 = abs(IKslow2 - max(IKslow2_trc));

    [~, tau_to_idx] = min(abs(Ito_hat*exp(-1) - Ito_trc));
    [~, tau1_idx] = min(abs(IKslow1_hat*exp(-1) - IKslow1_trc));
    [~, tau2_idx] = min(abs(IKslow2_hat*exp(-1) - IKslow2_trc));

    tau_to = abs(tau_to - t(tau_to_idx));
    tau_slow1 = abs(tau1 - t(tau1_idx));
    tau_slow2 = abs(tau2 - t(tau2_idx));

    z = trace_sim + amp_to + amp_slow1 + amp_slow2 + tau_to + tau_slow1 + tau_slow2;
end    
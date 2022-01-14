function z = Ktrace_fitness(X, ds_Ktrace, Iss)    
    [t, ~, A, ~] = Kv(X, -70, 450, 50, 25*1000, 50, 25*1000);

    % IKsum without Iss
    Ito_trc = A(:,5);
    IKslow1_trc = A(:,10);
    IKslow2_trc = A(:,15);
    IKsum = Ito_trc + IKslow1_trc + IKslow2_trc;

    % add the constant term (Iss) after the peak
    [~, peak_idx] = max(IKsum);
    IKsum(peak_idx:end) = IKsum(peak_idx:end) + Iss;

    trace_sim = dtw(IKsum, ds_Ktrace.I);
    disp(trace_sim);

    z = trace_sim;
end    
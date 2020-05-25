function z = Kv_anal_fitness(X, Ktrace, Iss)
    t = Ktrace.time;
    [~, peak_idx] = max(Ktrace.I);
    IKsum = Kv_anal(t(peak_idx:end), 50, X) + Iss;
    z = sum((IKsum - Ktrace.I(peak_idx:end)).^2);
end

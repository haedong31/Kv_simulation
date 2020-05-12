function e = IKslow2_fitness(X, Y, holding_p, holding_t, P1, P1_t, P2, P2_t)
    [~,~,A,~] = IKslow2(holding_p, holding_t, P1, P1_t, P2, P2_t, X);
    IKslow2_hat = A(:,65);
    e = sum((Y-max(IKslow2_hat)).^2);
end

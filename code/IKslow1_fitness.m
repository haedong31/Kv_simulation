function e = IKslow1_fitness(X, Y, holding_p, holding_t, P1, P1_t, P2, P2_t)
    [~,~,A,~] = IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, X);
    IKslow1_hat = A(:,67);
    e = sum((Y-IKslow1_hat(end)).^2);
end

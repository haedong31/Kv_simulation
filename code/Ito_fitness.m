function e = Ito_fitness(X, Y, holding_p, holding_t, P1, P1_t, P2, P2_t)
    [~,~,A,~] = Ito(holding_p, holding_t, P1, P1_t, P2, P2_t, X);
    Ito_hat = A(:,61);
    e = sum((Y-Ito_hat(end)).^2);
end
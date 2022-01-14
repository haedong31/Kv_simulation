function f = standard_exp(t, A, tau)
    f = A .* exp(-t ./ tau)';
end
function IKslow1 = IKslow1(t, V, X)
    Gto = 0.4067;
    Gslow = 0.16;
    ato0 = 0.265563e-2; % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977; % ito_f; Gating variable for transient outward K+ current
    aslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    islow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    
    % IKslow1
    aslow1_ss = 1./(1+exp(-(V-X(1))./X(2)));
    islow1_ss = 0.21./(1+exp(-(V-X(3))./X(4))) + 0.79./(1+exp(-(V-X(5))./X(6))); 
    tau_aslow1 = 0.4930 .* exp(-0.0629.*V) + X(7);
    tau_islow1 = 500.0 + X(8)./(1+exp((V+X(9))./0.0492));
    aslow1 = (aslow0 - aslow1_ss) .* exp(-t./tau_aslow1) + aslow1_ss;
    islow1 = (islow0 - islow1_ss) .* exp(-t./tau_islow1) + islow1_ss;
    IKslow1 = Gslow .* aslow1 .* islow1 .* (V+82.8);
end

function IKsum = Kv_anal(t, V, X)
    Gto = 0.4067;
    Gslow = 0.16;
    ato0 = 0.265563e-2; % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977; % ito_f; Gating variable for transient outward K+ current
    aslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    islow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    
    % Ito
    ato_ss = 1./(1+exp(-(V-X(1))./X(2)));
    ito_ss = 1./(1+exp((V-X(3))./X(4)));
    alpha_ato = 0.180640.*exp(0.0357700.*(V+30.0));
    beta_ato = 0.395600.*exp(-0.0623700.*(V+30.0));
    alpha_ito = (0.000152000.*exp((-V+13.5)./7.0))./(0.00670830.*exp((-V+33.5./7.0))+1.0);
    beta_ito = (0.000950000.*exp((V+X(5))./X(6)))./(0.0513350.*exp((V+X(5))./X(6))+1.0);
    tau_ato = 1./(alpha_ato + beta_ato);
    tau_ito = 1./(alpha_ito + beta_ito);
    ato = (ato0 - ato_ss) .* exp(-t./tau_ato) + ato_ss;
    ito = (ito0 - ito_ss) .* exp(-t./tau_ito) + ito_ss;
    Ito = Gto .* ato.^3 .* ito .* (V+82.8);
    
    % IKslow1
    aslow1_ss = 1./(1+exp(-(V-X(7))./X(8)));
    islow1_ss = 0.21./(1+exp(-(V-X(9))./X(10))) + 0.79./(1+exp(-(V-X(11))./X(12))); 
    tau_aslow1 = 0.4930 .* exp(-0.0629.*V) + X(15);
    tau_islow1 = 500.0 + X(13)./(1+exp((V+X(14))./0.0492));
    aslow1 = (aslow0 - aslow1_ss) .* exp(-t./tau_aslow1) + aslow1_ss;
    islow1 = (islow0 - islow1_ss) .* exp(-t./tau_islow1) + islow1_ss;
    IKslow1 = Gslow .* aslow1 .* islow1 .* (V+82.8);
    
    % IKslow2
    aslow2_ss = 1./(1+exp(-(V-X(16))./X(17)));
    islow2_ss = 0.21./(1+exp(-(V-X(18))./X(19))) + 0.79./(1+exp(-(V-X(20))./X(21))); 
    tau_aslow2 = 0.4930 .* exp(-0.0629.*V) + X(24);
    tau_islow2 = 500.0 + X(22)./(1+exp((V+X(23))./0.0492));
    aslow2 = (aslow0 - aslow2_ss) .* exp(-t./tau_aslow2) + aslow2_ss;
    islow2 = (islow0 - islow2_ss) .* exp(-t./tau_islow2) + islow2_ss;
    IKslow2 = Gslow .* aslow2 .* islow2 .* (V+82.8);
    
    IKsum = Ito + IKslow1 + IKslow2;
end

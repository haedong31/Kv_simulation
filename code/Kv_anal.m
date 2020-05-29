function IKsum = Kv_anal(t, V, X)
    Gto = 0.4067;
    Gslow = 0.16;
    ato0 = 0.265563e-2; % ato_f; Gating variable for transient outward K+ current
    ito0 = 0.999977; % ito_f; Gating variable for transient outward K+ current
    aslow0 = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    islow0 = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    
    % Ito; 20 control varialbes [1, 20]
    ato_ss = 1./(1+exp(-(V-X(1))./X(2)));
    ito_ss = 1./(1+exp((V-X(3))./X(4)));
    alpha_ato = X(5).*exp(X(6).*(V+X(7)));
    beta_ato = X(8).*exp(X(9).*(V+X(10)));
    alpha_ito = (X(11).*exp(-(V+X(12))./X(13)))./(X(14).*exp(-(V+X(15))./X(16))+1.0);
    beta_ito = (X(17).*exp((V+X(18))./X(19)))./(X(20).*exp((V+X(18))./X(19))+1.0);

    tau_ato = 1./(alpha_ato + beta_ato);
    tau_ito = 1./(alpha_ito + beta_ito);
    ato = (ato0 - ato_ss) .* exp(-t./tau_ato) + ato_ss;
    ito = (ito0 - ito_ss) .* exp(-t./tau_ito) + ito_ss;
    Ito = Gto .* ato.^3 .* ito .* (V+82.8);
    
    % IKslow1; 12 control varialbes [21, 32]
    aslow1_ss = 1./(1+exp(-(V+X(21))./X(22)));
    islow1_ss = 0.21./(1+exp(-(V+X(23))./X(24))) + 0.79./(1+exp(-(V+X(25))./X(26))); 
    tau_aslow1 = X(27) .* exp(X(28).*V) + X(29);
    tau_islow1 = 500.0 + X(30)./(1+exp((V+X(31))./X(32)));
    
    aslow1 = (aslow0 - aslow1_ss) .* exp(-t./tau_aslow1) + aslow1_ss;
    islow1 = (islow0 - islow1_ss) .* exp(-t./tau_islow1) + islow1_ss;
    IKslow1 = Gslow .* aslow1 .* islow1 .* (V+82.8);
    
    % IKslow2; 12 control variables [33, 44]
    aslow2_ss = 1./(1+exp(-(V+X(33))./X(34)));
    islow2_ss = 0.21./(1+exp(-(V+X(35))./X(36))) + 0.79./(1+exp(-(V+X(37))./X(38))); 
    tau_aslow2 = X(39) .* exp(X(40).*V) + X(41);
    tau_islow2 = 500.0 + X(42)./(1+exp((V+X(43))./X(44)));
    
    aslow2 = (aslow0 - aslow2_ss) .* exp(-t./tau_aslow2) + aslow2_ss;
    islow2 = (islow0 - islow2_ss) .* exp(-t./tau_islow2) + islow2_ss;
    IKslow2 = Gslow .* aslow2 .* islow2 .* (V+82.8);
    
    IKsum = Ito + IKslow1 + IKslow2;
end

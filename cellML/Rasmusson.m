
function [t, STATES, ALGEBRAIC, CONSTANTS] = Rasmusson()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % There are a total of 41 entries in each of the rate and state variable arrays.
    % There are a total of 73 entries in the constant variable array.
    algebraicVariableCount =72;
end

function [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 100];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [t, STATES] = ode15s(@(t, STATES)computeRates(t, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t);
end

function [STATES, CONSTANTS] = initConsts()
    CONSTANTS = []; STATES = [];

    % Table 8. Initial conditions
    STATES(:,1) = -82.4202;  % V; Membrane potential:mV
    STATES(:,2) = 0.115001;  % Cai; Myoplasmic Ca+ concentration:uM
    STATES(:,3) = 0.115001;  % Cass; Subspace SR Ca+ concentration:uM
    STATES(:,4) = 1299.5;  % CaJSR; JSR Ca+ concentration:uM
    STATES(:,5) = 1299.5;  % CaNSR; NSR Ca+ concentration:uM
    STATES(:,6) = 0;  % PRyR; RyR modulation factor
    STATES(:,7) = 11.2684;  % LTRPNCa; Concentration Ca+ bound low-affinity troponin-binding sites:uM
    STATES(:,8) = 125.29;  % HTEPNCa; Concentration Ca+ bound high-affinity troponin-binding sites:uM
    STATES(:,9) = 0.149102e-4;  % PO1; %Fraction of RyR channels in sate POl
    STATES(:,10) = 0.951726e-10;  % PO2; %Fraction of RyR channels in sate POl
    STATES(:,11) = 0.16774e-3;  % PC2; Fraction of RyR channels in sate PC2
    STATES(:,12) = 0.930308e-18;  % O; L-type Ca+ channel conducting state
    STATES(:,13) = 0.124216e-3;  % C2; L-type Ca+ channel closed state
    STATES(:,14) = 0.578679e-8;  % C3; L-type Ca+ channel closed state
    STATES(:,15) = 0.119816e-12;  % C4; L-type Ca+ channel closed state
    STATES(:,16) = 0.497923e-18;  % I1; L-type Ca+ channel inacticated state
    STATES(:,17) = 0.345847e-13;  % I2; L-type Ca+ channel inacticated state
    STATES(:,18) = 0.185106e-13;  % I3; L-type Ca+ channel inacticated state
    STATES(:,19) = 14237.1;  % Nai; Myoplasmic Na+ concentration:uM
    STATES(:,20) = 0.713483e-6;  % ONa; Open state of fast Na+ channel
    STATES(:,21) = 0.279132e-3;  % CNa1; Closed state of fast Na+ channel
    STATES(:,22) = 0.020752;  % CNa2; Closed state of fast Na+ channel
    STATES(:,23) = 0.673345e-6;  % I1Na; Slow inactivated state 1 of fast Na+ channel
    STATES(:,24) = 0.155787e-8;  % I2Na; Slow inactivated state 2 of fast Na+ channel
    STATES(:,25) = 0.153176e-3;  % IFNa; Fast inactivated state of fast Na+ channel
    STATES(:,26) = 0.0113879;  % ICNa2; Cloesd-inactivated state of fast Na+ channel
    STATES(:,27) = 0.34278;  % ICNa3; Cloesd-inactivated state of fast Na+ channel
    STATES(:,28) = 143720;  % Ki; Myoplasmic K+ concentration:uM
    STATES(:,29) = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    STATES(:,30) = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    STATES(:,31) = 0.417069e-3;  % ato_s; Gating variable for transient outward K+ current
    STATES(:,32) = 0.998543;  % ito_s; Gating variable for transient outward K+ current
    STATES(:,33) = 0.262753e-3;  % nKs; Gating variable for slow delayed-rectifier K+ current
    STATES(:,34) = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,35) = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,36) = 0.417069e-3;  % aKss; Gating variable for noninactivating steady-state K+ current
    STATES(:,37) = 1;  % iKss; Gating variable for noninactivating steady-state K+ current
    STATES(:,38) = 0.175298e-3;  % OK; mERG channel open state
    STATES(:,39) = 0.992513e-3;  % CK1; mERG channel closed state
    STATES(:,40) = 0.641229e-3;  % CK2; mERG channel closed state
    STATES(:,41) = 0.319129e-4;  % IK; mERG channel inactivated state
    
    % Table 2. Cell geometry parameters
    CONSTANTS(:,2) = 25.84e-6;  % Vmyo; Myoplasmic volume:ul
    CONSTANTS(:,3) = 0.12e-6;  % VJS; Juncional SR volume:ul
    CONSTANTS(:,4) = 2.098e-6;  % VNSR; Network SR volume:ul
    CONSTANTS(:,5) = 1.485e-9;  % Vss; Subspace volume:ul
    CONSTANTS(:,6) = 1.534e-4;  % Acap; Capacitive membrane area:cm^2
    
    % Table 3. Extracellular ion concentrations
    CONSTANTS(:,7) = 5400;  % Ko; Exracellular K+ concentration:uM
    CONSTANTS(:,8) = 140000;  % Nao; Exracellular Na+ concentration:uM
    CONSTANTS(:,9) = 1800;  % Cao; Exracellular Ca+ concentration:uM
    
    % Table 4. SR parameters
    CONSTANTS(:,26) = 4.5;  % nu1; Maximum RyR channel Ca+ permeability:ms^-1
    CONSTANTS(:,28) = 1.74e-5;  % nu2; Ca+ leak rate constant from the NSR:ms^-1
    CONSTANTS(:,27) = 20;  % tautr; Time constant for transfer from NSR to JSR:ms
    CONSTANTS(:,29) = 8;  % tauxfer; Time constant for transfer from subspace to myoplasm:ms
    CONSTANTS(:,30) = 0.45;  % nu3; SR Ca+-ATPase maximum pump rate:uM/ms
    CONSTANTS(:,31) = 0.5;  % Km_up; Half-saturation constant for SR Ca+-ATPase pump:uM
    CONSTANTS(:,35) = 0.006075;  % ka1; RyR Pc1-Po1 rate constant:uM^-4/ms
    CONSTANTS(:,36) = 0.07125;  % ka2; RyR Po1-Pc1 rate constant:ms^-1
    CONSTANTS(:,37) = 0.00405;  % kb1; RyR Po1-Po2 rate constant:uM^-3/ms
    CONSTANTS(:,38) = 0.965;  % kb2; RyR Po2-Po1 rate constant:ms^-1
    CONSTANTS(:,39) = 0.009;  % kc1; RyR Po1-Pc2 rate constant:ms^-1
    CONSTANTS(:,40) = 0.0008;  % kc2; RyR Pc2-Po1 rate constant:ms^-1
    CONSTANTS(:,41) = 3;  % m; RyR Ca+ cooperativity parameter Po1-Po2
    CONSTANTS(:,42) = 4;  % n; RyR Ca+ cooperativity parameter Pc1-Po1
    
    % Table 5. L-type Ca+ channel parameters
    CONSTANTS(:,44) = 0.1729;  % GCaL; Specific maximum conductivity for L-type Ca+ channel:mS/uF
    CONSTANTS(:,43) = 63;  % ECa; Reversal potential for L-type Ca+ channel:mV
    CONSTANTS(:,46) = 0.23324;  % Kpc_max; Maximum time constant for Ca+-induced inactivation:ms^-1
    CONSTANTS(:,47) = 20;  % Kpc_half; Half-saturation constant for Ca+-induced inactivation:uM
    CONSTANTS(:,45) = 0.0005;  % Kpcb; Voltage-insensitive rate constant for inactivation:ms^-1
    CONSTANTS(:,34) = 7;  % ICaL_max; Normalization constant for L-type Ca+ current:pA/pF 
    
    % Table 6. Buffering parameters
    CONSTANTS(:,32) = 70;  % LTRPN_tot; Total myoplasmic troponin low-affinity site concentration:uM
    CONSTANTS(:,33) = 140;  % HTRPN_tot; Total myoplasmic troponin high-affinity site concentration:uM
    CONSTANTS(:,22) = 0.00237;  % khtrpn1; Ca+ on rate constant for troponin high-affinity sites:uM^-1/ms
    CONSTANTS(:,23) = 3.2e-5;  % khtrpn2; Ca+ off rate constant for troponin high-affinity sites:ms^-1
    CONSTANTS(:,24) = 0.0327;  % kltrpn1; Ca+ on rate constant for troponin low-affinity sites:uM^-1/ms
    CONSTANTS(:,25) = 0.0196;  % kltrpn2; Ca+ off rate constant for troponin low-affinity sites:ms^-1
    CONSTANTS(:,18) = 50;  % CMDN_tot; Total myoplasmic calmodulin concentration:uM
    CONSTANTS(:,19) = 15000;  % CSQN_tot; Total junctional SR calsequestrin concentration:uM
    CONSTANTS(:,20) = 0.238;  % Km_CMDN; Ca+ half-saturaion constant for calmodulin:uM
    CONSTANTS(:,21) = 800;  % Km_CSQN; Ca+ half-saturaion constant for calsequestrin:uM
    
    % Table 7. Membrane current parameters
    CONSTANTS(:,12) = 96.5;  % F; Faraday constant:C/mmol
    CONSTANTS(:,11) = 298;  % T; Absolute temperature:K
    CONSTANTS(:,10) = 8.314;  % R; Ideal gas constant:J*mol^-1*K^-1
    CONSTANTS(:,50) = 292.8;  % kNaCa; Scaling factor of Na+/Ca+ exchange:pA/pF
    CONSTANTS(:,51) = 87500;  % Km_Na; Na+ half-saturation constant for Na+/Ca+ exchange:uM
    CONSTANTS(:,52) = 1380;  %  Km_Ca; Ca+ half-saturation constant for Na+/Ca+ exchange:uM
    CONSTANTS(:,53) = 0.1;  % ksat; Na+/Ca+ exchange saaturation factor at very negative potentials
    CONSTANTS(:,54) = 0.35;  % yita; Contals voltage dependence of Na+/Ca+ exchange
    CONSTANTS(:,66) = 0.88;  % INaK_max; Maximum Na+/K+ exchange current:pA/pF 
    CONSTANTS(:,67) = 21000;  % Km_Nai; Na+ half-saturation constant for Na+/K+ exchange currant:uM
    CONSTANTS(:,68) = 1500;  % Km_Ko; K+ half-saturation constant for Na+/K+ exchange currant:uM
    CONSTANTS(:,48) = 1;  % Ip_Ca_max; Maximun Ca+ pump current:pA/pF
    CONSTANTS(:,49) = 0.5;  % Km_p_Ca; Ca+ half-saturation constant for Ca+pump current:uM
    CONSTANTS(:,55) = 0.000367;  % GCab; Maximun background Ca+ current conductance:mS/uF
    CONSTANTS(:,56) = 13;  % GNa; Maximun fast Na+ current conductance:mS/uF
    CONSTANTS(:,57) = 0.0026;  % GNab; Maximun background Na+ current conductance:mS/uF
    CONSTANTS(:,60) = 0.00575;  % GKs; Maximum slow delayed-rectifier K+ current conductance:mS/uF
    
    % Table 7. Parameters for apex
    CONSTANTS(:,58) = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    CONSTANTS(:,59) = 0;  % GKtos; Maximum transient outward K+ current conductance(apex):mS/uF
    CONSTANTS(:,61) = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    CONSTANTS(:,62) = 0.05;  % GKss; Maximum noninactivating steady-state K+ current conductance(apex):mS/uF
    
    % Table 7. Membrane current parameters
    CONSTANTS(:,63) = 0.078;  % GKr; Maximum rapid delayed-rectifier K+ current conductance:mS/uF
    CONSTANTS(:,65) = 0.023761;  % kf; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    CONSTANTS(:,64) = 0.036778;  % kb; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    CONSTANTS(:,69) = 10;  % GCl_Ca; Maximum Ca+-activated Cl- current conductance:mS/uF
    CONSTANTS(:,71) = 10;  % Km_Cl; Half-saturaon constant for Ca+-activated Cl- current:uM
    CONSTANTS(:,70) = -40;  % ECl; Reversal potential for Ca+-activated Cl- current:mV
    
    % Parameters not in the Simulink model
    CONSTANTS(:,1) = 1;
    CONSTANTS(:,13) = 20;
    CONSTANTS(:,14) = 100000;
    CONSTANTS(:,15) = 71.43;
    CONSTANTS(:,16) = 0.5;
    CONSTANTS(:,17) = -80;
    CONSTANTS(:,72) =  (1.00000./7.00000).*(exp(CONSTANTS(:,8)./67300.0) - 1.00000);
    CONSTANTS(:,72) = 0.00000;

    if (isempty(STATES)), warning('Initial values for states not set'); end
end

function [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    
    % externally applied voltage (voltage clamp)
    ALGEBRAIC(:,72) = arrayfun(@volt_clamp, t);
    
    % A94; iKss if CONSTANTS(:,72)=0  or A110; sigma if CONSTANTS(:,72)=A110
    RATES(:,37) = CONSTANTS(:,72); 
    % A16; LTRPNCa
    RATES(:,7) =  CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7)) -  CONSTANTS(:,25).*STATES(:,7);
    % A17; HTEPNCa
    RATES(:,8) =  CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8)) -  CONSTANTS(:,23).*STATES(:,8);
    % A20; PO2
    RATES(:,10) =  CONSTANTS(:,37).*power(STATES(:,3), CONSTANTS(:,41)).*STATES(:,9) -  CONSTANTS(:,38).*STATES(:,10);
    % A21; PC2
    RATES(:,11) =  CONSTANTS(:,39).*STATES(:,9) -  CONSTANTS(:,40).*STATES(:,11);
    % A19; PC1
    ALGEBRAIC(:,2) = 1.00000 - (STATES(:,11)+STATES(:,9)+STATES(:,10));
    % A18; PO1
    RATES(:,9) = ( CONSTANTS(:,35).*power(STATES(:,3), CONSTANTS(:,42)).*ALGEBRAIC(:,2)+ CONSTANTS(:,38).*STATES(:,10)+ CONSTANTS(:,40).*STATES(:,11)) - ( CONSTANTS(:,36).*STATES(:,9)+ CONSTANTS(:,37).*power(STATES(:,3), CONSTANTS(:,41)).*STATES(:,9)+ CONSTANTS(:,39).*STATES(:,9));
    % A71; alpha_a
    ALGEBRAIC(:,5) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,72)+30.0000));
    % A72; beta_a
    ALGEBRAIC(:,15) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,72)+30.0000));
    % A69; ato_f
    RATES(:,29) =  ALGEBRAIC(:,5).*(1.00000 - STATES(:,29)) -  ALGEBRAIC(:,15).*STATES(:,29);
    % A73; alpha_i
    ALGEBRAIC(:,6) = ( 0.000152000.*exp( - (ALGEBRAIC(:,72)+13.5000)./7.00000))./( 0.00670830.*exp( - (ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    % A74; beta_i
    ALGEBRAIC(:,16) = ( 0.000950000.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000))./( 0.0513350.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    % A70; ito_f
    RATES(:,30) =  ALGEBRAIC(:,6).*(1.00000 - STATES(:,30)) -  ALGEBRAIC(:,16).*STATES(:,30);
    % A78; a_ss
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,72)+22.5000)./7.70000));
    % A80; tau_tas
    ALGEBRAIC(:,17) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    % A76; ato_s
    RATES(:,31) = (ALGEBRAIC(:,7) - STATES(:,31))./ALGEBRAIC(:,17);
    % A79; i_ss
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A81; tau_tis
    ALGEBRAIC(:,18) = 270.000+1050.00./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A77; ito_s
    RATES(:,32) = (ALGEBRAIC(:,8) - STATES(:,32))./ALGEBRAIC(:,18);
    % A85; alpha_n
    ALGEBRAIC(:,9) = ( 4.81333e-06.*(ALGEBRAIC(:,72)+26.5000))./(1.00000 - exp(  - 0.128000.*(ALGEBRAIC(:,72)+26.5000)));
    % A86; beta_n
    ALGEBRAIC(:,19) =  9.53333e-05.*exp(  - 0.0380000.*(ALGEBRAIC(:,72)+26.5000));
    % A84; nKs
    RATES(:,33) =  ALGEBRAIC(:,9).*(1.00000 - STATES(:,33)) -  ALGEBRAIC(:,19).*STATES(:,33);
    % A90; tau_aur
    ALGEBRAIC(:,20) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    % A88; aur
    RATES(:,34) = (ALGEBRAIC(:,7) - STATES(:,34))./ALGEBRAIC(:,20);
    % A91; tau_iur
    ALGEBRAIC(:,21) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A89; iur
    RATES(:,35) = (ALGEBRAIC(:,8) - STATES(:,35))./ALGEBRAIC(:,21);
    % A95; tau_Kss
    ALGEBRAIC(:,22) =  39.3000.*exp(  - 0.0862000.*ALGEBRAIC(:,72))+13.1700;
    % A93; aKss
    RATES(:,36) = (ALGEBRAIC(:,7) - STATES(:,36))./ALGEBRAIC(:,22);
    % A104; alpha_a1
    ALGEBRAIC(:,11) =  0.0137330.*exp( 0.0381980.*ALGEBRAIC(:,72));
    % A105; beta_a1
    ALGEBRAIC(:,24) =  6.89000e-05.*exp(  - 0.0417800.*ALGEBRAIC(:,72));
    % A99; CK2
    RATES(:,40) = ( CONSTANTS(:,65).*STATES(:,39)+ ALGEBRAIC(:,24).*STATES(:,38)) - ( CONSTANTS(:,64).*STATES(:,40)+ ALGEBRAIC(:,11).*STATES(:,40));
    % A24; C1
    ALGEBRAIC(:,3) = 1.00000 - (STATES(:,12)+STATES(:,13)+STATES(:,14)+STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,18));
    % A31; alpha
    ALGEBRAIC(:,13) = ( 0.400000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000).*((1.00000+ 0.700000.*exp( - power(ALGEBRAIC(:,72)+40.0000, 2.00000)./10.0000)) -  0.750000.*exp( - power(ALGEBRAIC(:,72)+20.0000, 2.00000)./400.000)))./(1.00000+ 0.120000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000));
    % A32; beta
    ALGEBRAIC(:,26) =  0.0500000.*exp( - (ALGEBRAIC(:,72)+12.0000)./13.0000);
    % A25; C2
    RATES(:,13) = ( 4.00000.*ALGEBRAIC(:,13).*ALGEBRAIC(:,3)+ 2.00000.*ALGEBRAIC(:,26).*STATES(:,14)) - ( ALGEBRAIC(:,26).*STATES(:,13)+ 3.00000.*ALGEBRAIC(:,13).*STATES(:,13));
    % A26; C3
    RATES(:,14) = ( 3.00000.*ALGEBRAIC(:,13).*STATES(:,13)+ 3.00000.*ALGEBRAIC(:,26).*STATES(:,15)) - ( 2.00000.*ALGEBRAIC(:,26).*STATES(:,14)+ 2.00000.*ALGEBRAIC(:,13).*STATES(:,14));
    % A97; CK0
    ALGEBRAIC(:,10) = 1.00000 - (STATES(:,39)+STATES(:,40)+STATES(:,38)+STATES(:,41));
    % A102; alpha_a0
    ALGEBRAIC(:,23) =  0.0223480.*exp( 0.0117600.*ALGEBRAIC(:,72));
    % A103; beta_a0
    ALGEBRAIC(:,28) =  0.0470020.*exp(  - 0.0631000.*ALGEBRAIC(:,72));
    % A98; CK1
    RATES(:,39) = ( ALGEBRAIC(:,23).*ALGEBRAIC(:,10)+ CONSTANTS(:,64).*STATES(:,40)) - ( ALGEBRAIC(:,28).*STATES(:,39)+ CONSTANTS(:,65).*STATES(:,39));
    % A106; alpha_i
    ALGEBRAIC(:,29) =  0.0908210.*exp( 0.0233910.*(ALGEBRAIC(:,72)+5.00000));
    % A107; beta_i
    ALGEBRAIC(:,33) =  0.00649700.*exp(  - 0.0326800.*(ALGEBRAIC(:,72)+5.00000));
    % A100; OK
    RATES(:,38) = ( ALGEBRAIC(:,11).*STATES(:,40)+ ALGEBRAIC(:,33).*STATES(:,41)) - ( ALGEBRAIC(:,24).*STATES(:,38)+ ALGEBRAIC(:,29).*STATES(:,38));
    % A101; IK
    RATES(:,41) =  ALGEBRAIC(:,29).*STATES(:,38) -  ALGEBRAIC(:,33).*STATES(:,41);
    % A33; gamma
    ALGEBRAIC(:,31) = ( CONSTANTS(:,46).*STATES(:,3))./(CONSTANTS(:,47)+STATES(:,3));
    % A34; K_pcf
    ALGEBRAIC(:,35) =  13.0000.*(1.00000 - exp( - power(ALGEBRAIC(:,72)+14.5000, 2.00000)./100.000));
    % A23; O
    RATES(:,12) = ( ALGEBRAIC(:,13).*STATES(:,15)+ CONSTANTS(:,45).*STATES(:,16)+ 0.00100000.*( ALGEBRAIC(:,13).*STATES(:,17) -  ALGEBRAIC(:,35).*STATES(:,12))) - ( 4.00000.*ALGEBRAIC(:,26).*STATES(:,12)+ ALGEBRAIC(:,31).*STATES(:,12));
    % A27; C4
    RATES(:,15) = ( 2.00000.*ALGEBRAIC(:,13).*STATES(:,14)+ 4.00000.*ALGEBRAIC(:,26).*STATES(:,12)+ 0.0100000.*( 4.00000.*CONSTANTS(:,45).*ALGEBRAIC(:,26).*STATES(:,16) -  ALGEBRAIC(:,13).*ALGEBRAIC(:,31).*STATES(:,15))+ 0.00200000.*( 4.00000.*ALGEBRAIC(:,26).*STATES(:,17) -  ALGEBRAIC(:,35).*STATES(:,15))+ 4.00000.*ALGEBRAIC(:,26).*CONSTANTS(:,45).*STATES(:,18)) - ( 3.00000.*ALGEBRAIC(:,26).*STATES(:,15)+ ALGEBRAIC(:,13).*STATES(:,15)+ 1.00000.*ALGEBRAIC(:,31).*ALGEBRAIC(:,35).*STATES(:,15));
    % A28; I1
    RATES(:,16) = ( ALGEBRAIC(:,31).*STATES(:,12)+ 0.00100000.*( ALGEBRAIC(:,13).*STATES(:,18) -  ALGEBRAIC(:,35).*STATES(:,16))+ 0.0100000.*( ALGEBRAIC(:,13).*ALGEBRAIC(:,31).*STATES(:,15) -  4.00000.*ALGEBRAIC(:,26).*ALGEBRAIC(:,35).*STATES(:,16))) -  CONSTANTS(:,45).*STATES(:,16);
    % A29; I2
    RATES(:,17) = ( 0.00100000.*( ALGEBRAIC(:,35).*STATES(:,12) -  ALGEBRAIC(:,13).*STATES(:,17))+ CONSTANTS(:,45).*STATES(:,18)+ 0.00200000.*( ALGEBRAIC(:,35).*STATES(:,15) -  4.00000.*ALGEBRAIC(:,26).*STATES(:,17))) -  ALGEBRAIC(:,31).*STATES(:,17);
    % A30; I3
    RATES(:,18) = ( 0.00100000.*( ALGEBRAIC(:,35).*STATES(:,16) -  ALGEBRAIC(:,13).*STATES(:,18))+ ALGEBRAIC(:,31).*STATES(:,17)+ 1.00000.*ALGEBRAIC(:,31).*ALGEBRAIC(:,35).*STATES(:,15)) - ( 4.00000.*ALGEBRAIC(:,26).*CONSTANTS(:,45).*STATES(:,18)+ CONSTANTS(:,45).*STATES(:,18));
    % A8; B_JSR
    ALGEBRAIC(:,30) = power(1.00000+( CONSTANTS(:,19).*CONSTANTS(:,21))./power(CONSTANTS(:,21)+STATES(:,4), 2.00000),  - 1.00000);
    % A9; J_rel
    ALGEBRAIC(:,34) =  CONSTANTS(:,26).*(STATES(:,9)+STATES(:,10)).*(STATES(:,4) - STATES(:,3)).*STATES(:,6);
    % A10; J_tr
    ALGEBRAIC(:,37) = (STATES(:,5) - STATES(:,4))./CONSTANTS(:,27);
    % A4; CaJSR
    RATES(:,4) =  ALGEBRAIC(:,30).*(ALGEBRAIC(:,37) - ALGEBRAIC(:,34));
    % A12; J_leak
    ALGEBRAIC(:,41) =  CONSTANTS(:,28).*(STATES(:,5) - STATES(:,2));
    % A13; J_up
    ALGEBRAIC(:,43) = ( CONSTANTS(:,30).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,31), 2.00000)+power(STATES(:,2), 2.00000));
    % A5; CaNSR
    RATES(:,5) = ( (ALGEBRAIC(:,43) - ALGEBRAIC(:,41)).*CONSTANTS(:,2))./CONSTANTS(:,4) - ( ALGEBRAIC(:,37).*CONSTANTS(:,3))./CONSTANTS(:,4);
    % A42; CNa3
    ALGEBRAIC(:,4) = 1.00000 - (STATES(:,20)+STATES(:,21)+STATES(:,22)+STATES(:,25)+STATES(:,23)+STATES(:,24)+STATES(:,26)+STATES(:,27));
    % A51; alpha_Na11
    ALGEBRAIC(:,14) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./17.0000)+ 0.200000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A54; beta_Na11
    ALGEBRAIC(:,36) =  0.191700.*exp( - (ALGEBRAIC(:,72)+2.50000)./20.3000);
    % A52; alpha_Na12
    ALGEBRAIC(:,27) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./15.0000)+ 0.230000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A55; beta_Na12
    ALGEBRAIC(:,38) =  0.200000.*exp( - (ALGEBRAIC(:,72) - 2.50000)./20.3000);
    % A57; alpha_Na3
    ALGEBRAIC(:,42) =  7.00000e-07.*exp( - (ALGEBRAIC(:,72)+7.00000)./7.70000);
    % A58; beta_Na3
    ALGEBRAIC(:,44) = 0.00840000+ 2.00000e-05.*(ALGEBRAIC(:,72)+7.00000);
    % A43; CNa2
    RATES(:,22) = ( ALGEBRAIC(:,14).*ALGEBRAIC(:,4)+ ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,42).*STATES(:,26)) - ( ALGEBRAIC(:,36).*STATES(:,22)+ ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,44).*STATES(:,22));
    % A53; alpha_Na13
    ALGEBRAIC(:,32) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./12.0000)+ 0.250000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A56; beta_Na13
    ALGEBRAIC(:,40) =  0.220000.*exp( - (ALGEBRAIC(:,72) - 7.50000)./20.3000);
    % A44; CNa1
    RATES(:,21) = ( ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,42).*STATES(:,25)) - ( ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,44).*STATES(:,21));
    % A49; ICNa2
    RATES(:,26) = ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,38).*STATES(:,25)+ ALGEBRAIC(:,44).*STATES(:,22)) - ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,27).*STATES(:,26)+ ALGEBRAIC(:,42).*STATES(:,26));
    % A50; ICNa3
    RATES(:,27) = ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,44).*ALGEBRAIC(:,4)) - ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,42).*STATES(:,27));
    % A22; I_CaL
    ALGEBRAIC(:,47) =  CONSTANTS(:,44).*STATES(:,12).*(ALGEBRAIC(:,72) - CONSTANTS(:,43));
    % A7; B_ss
    ALGEBRAIC(:,25) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,3), 2.00000),  - 1.00000);
    % A11; J_xfer
    ALGEBRAIC(:,39) = (STATES(:,3) - STATES(:,2))./CONSTANTS(:,29);
    % A3; Cass
    RATES(:,3) =  ALGEBRAIC(:,25).*(( ALGEBRAIC(:,34).*CONSTANTS(:,3))./CONSTANTS(:,5) - (( ALGEBRAIC(:,39).*CONSTANTS(:,2))./CONSTANTS(:,5)+( ALGEBRAIC(:,47).*CONSTANTS(:,6).*CONSTANTS(:,1))./( 2.00000.*CONSTANTS(:,5).*CONSTANTS(:,12))));
    % A15; PRyR
    RATES(:,6) =   - 0.0400000.*STATES(:,6) -  (( 0.100000.*ALGEBRAIC(:,47))./CONSTANTS(:,34)).*exp( - power(ALGEBRAIC(:,72) - 5.00000, 2.00000)./648.000);
    % A59; alpha_Na2
    ALGEBRAIC(:,46) = 1.00000./( 0.188495.*exp( - (ALGEBRAIC(:,72)+7.00000)./16.6000)+0.393956);
    % A60; beta_Na2
    ALGEBRAIC(:,48) = ( ALGEBRAIC(:,32).*ALGEBRAIC(:,46).*ALGEBRAIC(:,42))./( ALGEBRAIC(:,40).*ALGEBRAIC(:,44));
    % A45; ONa
    RATES(:,20) = ( ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,48).*STATES(:,25)) - ( ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,46).*STATES(:,20));
    % A61; alpha_Na4
    ALGEBRAIC(:,50) = ALGEBRAIC(:,46)./1000.00;
    % A62; beta_Na4
    ALGEBRAIC(:,52) = ALGEBRAIC(:,42);
    % A46; IFNa
    RATES(:,25) = ( ALGEBRAIC(:,46).*STATES(:,20)+ ALGEBRAIC(:,44).*STATES(:,21)+ ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,27).*STATES(:,26)) - ( ALGEBRAIC(:,48).*STATES(:,25)+ ALGEBRAIC(:,42).*STATES(:,25)+ ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,38).*STATES(:,25));
    % A35; I_p(Ca)
    ALGEBRAIC(:,49) = ( CONSTANTS(:,48).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,49), 2.00000)+power(STATES(:,2), 2.00000));
    % A36; I_NaCa
    ALGEBRAIC(:,51) =  (( (( (( CONSTANTS(:,50).*1.00000)./(power(CONSTANTS(:,51), 3.00000)+power(CONSTANTS(:,8), 3.00000))).*1.00000)./(CONSTANTS(:,52)+CONSTANTS(:,9))).*1.00000)./(1.00000+ CONSTANTS(:,53).*exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*( exp(( CONSTANTS(:,54).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,19), 3.00000).*CONSTANTS(:,9) -  exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8), 3.00000).*STATES(:,2));
    % A38; E_CaN
    ALGEBRAIC(:,53) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./( 2.00000.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./STATES(:,2));
    % A37; I_Cab
    ALGEBRAIC(:,55) =  CONSTANTS(:,55).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,53));
    % A6; Bi
    ALGEBRAIC(:,12) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,2), 2.00000),  - 1.00000);
    % A14; J_trpn
    ALGEBRAIC(:,45) = ( CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8))+ CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7))) - ( CONSTANTS(:,23).*STATES(:,8)+ CONSTANTS(:,25).*STATES(:,7));
    % A2; Cai
    RATES(:,2) =  ALGEBRAIC(:,12).*((ALGEBRAIC(:,41)+ALGEBRAIC(:,39)) - (ALGEBRAIC(:,43)+ALGEBRAIC(:,45)+( ((ALGEBRAIC(:,55)+ALGEBRAIC(:,49)) -  2.00000.*ALGEBRAIC(:,51)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( 2.00000.*CONSTANTS(:,2).*CONSTANTS(:,12))));
    % A63; alpha_Na5
    ALGEBRAIC(:,54) = ALGEBRAIC(:,46)./95000.0;
    % A64; beta_Na5
    ALGEBRAIC(:,56) = ALGEBRAIC(:,42)./50.0000;
    % A47; I1Na
    RATES(:,23) = ( ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,56).*STATES(:,24)) - ( ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,54).*STATES(:,23));
    % A48; I2Na
    RATES(:,24) =  ALGEBRAIC(:,54).*STATES(:,23) -  ALGEBRAIC(:,56).*STATES(:,24);
    % A41; E_Na
    ALGEBRAIC(:,57) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.900000.*CONSTANTS(:,8)+ 0.100000.*CONSTANTS(:,7))./( 0.900000.*STATES(:,19)+ 0.100000.*STATES(:,28)));
    % A40; I_Na
    ALGEBRAIC(:,58) =  CONSTANTS(:,56).*STATES(:,20).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    % A65; I_Nab
    ALGEBRAIC(:,59) =  CONSTANTS(:,57).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    % A109; f_NaK
    ALGEBRAIC(:,68) = 1.00000./(1.00000+ 0.124500.*exp((  - 0.100000.*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11)))+ 0.0365000.*CONSTANTS(:,72).*exp((  - ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))));
    % A108; I_NaK
    ALGEBRAIC(:,69) = ( (( CONSTANTS(:,66).*ALGEBRAIC(:,68).*1.00000)./(1.00000+power(CONSTANTS(:,67)./STATES(:,19), 1.50000))).*CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,68));
    % A39; Nai
    RATES(:,19) = (  - (ALGEBRAIC(:,58)+ALGEBRAIC(:,59)+ 3.00000.*ALGEBRAIC(:,69)+ 3.00000.*ALGEBRAIC(:,51)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( CONSTANTS(:,2).*CONSTANTS(:,12));
    % A68; E_K
    ALGEBRAIC(:,60) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,28));
    % A67; I_Kto,f  
    ALGEBRAIC(:,61) =  CONSTANTS(:,58).*power(STATES(:,29), 3.00000).*STATES(:,30).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A75; I_Kto,s
    ALGEBRAIC(:,62) =  CONSTANTS(:,59).*STATES(:,31).*STATES(:,32).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A82; I_K1
    ALGEBRAIC(:,63) = ( (( 0.293800.*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.000)).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60)))./(1.00000+exp( 0.0896000.*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60))));
    % A83; I_Ks
    ALGEBRAIC(:,64) =  CONSTANTS(:,60).*power(STATES(:,33), 2.00000).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % 87; I_kUR
    ALGEBRAIC(:,65) =  CONSTANTS(:,61).*STATES(:,34).*STATES(:,35).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A92; I_Kss
    ALGEBRAIC(:,66) =  CONSTANTS(:,62).*STATES(:,36).*STATES(:,37).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A96; I_Kr
    ALGEBRAIC(:,67) =  CONSTANTS(:,63).*STATES(:,38).*(ALGEBRAIC(:,72) -  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.980000.*CONSTANTS(:,7)+ 0.0200000.*CONSTANTS(:,8))./( 0.980000.*STATES(:,28)+ 0.0200000.*STATES(:,19))));
    % A66; Ki
    RATES(:,28) = (  - ((ALGEBRAIC(:,61)+ALGEBRAIC(:,62)+ALGEBRAIC(:,63)+ALGEBRAIC(:,64)+ALGEBRAIC(:,66)+ALGEBRAIC(:,65)+ALGEBRAIC(:,67)) -  2.00000.*ALGEBRAIC(:,69)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( CONSTANTS(:,2).*CONSTANTS(:,12));
    % I_stim
    ALGEBRAIC(:,1) = piecewise({t>=CONSTANTS(:,13)&t<=CONSTANTS(:,14)&(t - CONSTANTS(:,13)) -  floor((t - CONSTANTS(:,13))./CONSTANTS(:,15)).*CONSTANTS(:,15)<=CONSTANTS(:,16), CONSTANTS(:,17) }, 0.00000);
    % A112; O_ClCa
    ALGEBRAIC(:,70) = 0.200000./(1.00000+exp( - (ALGEBRAIC(:,72) - 46.7000)./7.80000));
    % A111; i_ClCa
    ALGEBRAIC(:,71) =  (( CONSTANTS(:,69).*ALGEBRAIC(:,70).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,71))).*(ALGEBRAIC(:,72) - CONSTANTS(:,70));
    % A1; V
    RATES(:,1) =  0;
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t)
    ALGEBRAIC(:,2) = 1.00000 - (STATES(:,11)+STATES(:,9)+STATES(:,10));
    ALGEBRAIC(:,5) =  0.180640.*exp( 0.0357700.*(STATES(:,1)+30.0000));
    ALGEBRAIC(:,15) =  0.395600.*exp(  - 0.0623700.*(STATES(:,1)+30.0000));
    ALGEBRAIC(:,6) = ( 0.000152000.*exp( - (STATES(:,1)+13.5000)./7.00000))./( 0.00670830.*exp( - (STATES(:,1)+33.5000)./7.00000)+1.00000);
    ALGEBRAIC(:,16) = ( 0.000950000.*exp((STATES(:,1)+33.5000)./7.00000))./( 0.0513350.*exp((STATES(:,1)+33.5000)./7.00000)+1.00000);
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp( - (STATES(:,1)+22.5000)./7.70000));
    ALGEBRAIC(:,17) =  0.493000.*exp(  - 0.0629000.*STATES(:,1))+2.05800;
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+45.2000)./5.70000));
    ALGEBRAIC(:,18) = 270.000+1050.00./(1.00000+exp((STATES(:,1)+45.2000)./5.70000));
    ALGEBRAIC(:,9) = ( 4.81333e-06.*(STATES(:,1)+26.5000))./(1.00000 - exp(  - 0.128000.*(STATES(:,1)+26.5000)));
    ALGEBRAIC(:,19) =  9.53333e-05.*exp(  - 0.0380000.*(STATES(:,1)+26.5000));
    ALGEBRAIC(:,20) =  0.493000.*exp(  - 0.0629000.*STATES(:,1))+2.05800;
    ALGEBRAIC(:,21) = 1200.00 - 170.000./(1.00000+exp((STATES(:,1)+45.2000)./5.70000));
    ALGEBRAIC(:,22) =  39.3000.*exp(  - 0.0862000.*STATES(:,1))+13.1700;
    ALGEBRAIC(:,11) =  0.0137330.*exp( 0.0381980.*STATES(:,1));
    ALGEBRAIC(:,24) =  6.89000e-05.*exp(  - 0.0417800.*STATES(:,1));
    ALGEBRAIC(:,3) = 1.00000 - (STATES(:,12)+STATES(:,13)+STATES(:,14)+STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,18));
    ALGEBRAIC(:,13) = ( 0.400000.*exp((STATES(:,1)+12.0000)./10.0000).*((1.00000+ 0.700000.*exp( - power(STATES(:,1)+40.0000, 2.00000)./10.0000)) -  0.750000.*exp( - power(STATES(:,1)+20.0000, 2.00000)./400.000)))./(1.00000+ 0.120000.*exp((STATES(:,1)+12.0000)./10.0000));
    ALGEBRAIC(:,26) =  0.0500000.*exp( - (STATES(:,1)+12.0000)./13.0000);
    ALGEBRAIC(:,10) = 1.00000 - (STATES(:,39)+STATES(:,40)+STATES(:,38)+STATES(:,41));
    ALGEBRAIC(:,23) =  0.0223480.*exp( 0.0117600.*STATES(:,1));
    ALGEBRAIC(:,28) =  0.0470020.*exp(  - 0.0631000.*STATES(:,1));
    ALGEBRAIC(:,29) =  0.0908210.*exp( 0.0233910.*(STATES(:,1)+5.00000));
    ALGEBRAIC(:,33) =  0.00649700.*exp(  - 0.0326800.*(STATES(:,1)+5.00000));
    ALGEBRAIC(:,31) = ( CONSTANTS(:,46).*STATES(:,3))./(CONSTANTS(:,47)+STATES(:,3));
    ALGEBRAIC(:,35) =  13.0000.*(1.00000 - exp( - power(STATES(:,1)+14.5000, 2.00000)./100.000));
    ALGEBRAIC(:,30) = power(1.00000+( CONSTANTS(:,19).*CONSTANTS(:,21))./power(CONSTANTS(:,21)+STATES(:,4), 2.00000),  - 1.00000);
    ALGEBRAIC(:,34) =  CONSTANTS(:,26).*(STATES(:,9)+STATES(:,10)).*(STATES(:,4) - STATES(:,3)).*STATES(:,6);
    ALGEBRAIC(:,37) = (STATES(:,5) - STATES(:,4))./CONSTANTS(:,27);
    ALGEBRAIC(:,41) =  CONSTANTS(:,28).*(STATES(:,5) - STATES(:,2));
    ALGEBRAIC(:,43) = ( CONSTANTS(:,30).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,31), 2.00000)+power(STATES(:,2), 2.00000));
    ALGEBRAIC(:,4) = 1.00000 - (STATES(:,20)+STATES(:,21)+STATES(:,22)+STATES(:,25)+STATES(:,23)+STATES(:,24)+STATES(:,26)+STATES(:,27));
    ALGEBRAIC(:,14) = 3.80200./( 0.102700.*exp( - (STATES(:,1)+2.50000)./17.0000)+ 0.200000.*exp( - (STATES(:,1)+2.50000)./150.000));
    ALGEBRAIC(:,36) =  0.191700.*exp( - (STATES(:,1)+2.50000)./20.3000);
    ALGEBRAIC(:,27) = 3.80200./( 0.102700.*exp( - (STATES(:,1)+2.50000)./15.0000)+ 0.230000.*exp( - (STATES(:,1)+2.50000)./150.000));
    ALGEBRAIC(:,38) =  0.200000.*exp( - (STATES(:,1) - 2.50000)./20.3000);
    ALGEBRAIC(:,42) =  7.00000e-07.*exp( - (STATES(:,1)+7.00000)./7.70000);
    ALGEBRAIC(:,44) = 0.00840000+ 2.00000e-05.*(STATES(:,1)+7.00000);
    ALGEBRAIC(:,32) = 3.80200./( 0.102700.*exp( - (STATES(:,1)+2.50000)./12.0000)+ 0.250000.*exp( - (STATES(:,1)+2.50000)./150.000));
    ALGEBRAIC(:,40) =  0.220000.*exp( - (STATES(:,1) - 7.50000)./20.3000);
    ALGEBRAIC(:,47) =  CONSTANTS(:,44).*STATES(:,12).*(STATES(:,1) - CONSTANTS(:,43));
    ALGEBRAIC(:,25) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,3), 2.00000),  - 1.00000);
    ALGEBRAIC(:,39) = (STATES(:,3) - STATES(:,2))./CONSTANTS(:,29);
    ALGEBRAIC(:,46) = 1.00000./( 0.188495.*exp( - (STATES(:,1)+7.00000)./16.6000)+0.393956);
    ALGEBRAIC(:,48) = ( ALGEBRAIC(:,32).*ALGEBRAIC(:,46).*ALGEBRAIC(:,42))./( ALGEBRAIC(:,40).*ALGEBRAIC(:,44));
    ALGEBRAIC(:,50) = ALGEBRAIC(:,46)./1000.00;
    ALGEBRAIC(:,52) = ALGEBRAIC(:,42);
    ALGEBRAIC(:,49) = ( CONSTANTS(:,48).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,49), 2.00000)+power(STATES(:,2), 2.00000));
    ALGEBRAIC(:,51) =  (( (( (( CONSTANTS(:,50).*1.00000)./(power(CONSTANTS(:,51), 3.00000)+power(CONSTANTS(:,8), 3.00000))).*1.00000)./(CONSTANTS(:,52)+CONSTANTS(:,9))).*1.00000)./(1.00000+ CONSTANTS(:,53).*exp(( (CONSTANTS(:,54) - 1.00000).*STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*( exp(( CONSTANTS(:,54).*STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,19), 3.00000).*CONSTANTS(:,9) -  exp(( (CONSTANTS(:,54) - 1.00000).*STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8), 3.00000).*STATES(:,2));
    ALGEBRAIC(:,53) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./( 2.00000.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./STATES(:,2));
    ALGEBRAIC(:,55) =  CONSTANTS(:,55).*(STATES(:,1) - ALGEBRAIC(:,53));
    ALGEBRAIC(:,12) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,2), 2.00000),  - 1.00000);
    ALGEBRAIC(:,45) = ( CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8))+ CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7))) - ( CONSTANTS(:,23).*STATES(:,8)+ CONSTANTS(:,25).*STATES(:,7));
    ALGEBRAIC(:,54) = ALGEBRAIC(:,46)./95000.0;
    ALGEBRAIC(:,56) = ALGEBRAIC(:,42)./50.0000;
    ALGEBRAIC(:,57) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.900000.*CONSTANTS(:,8)+ 0.100000.*CONSTANTS(:,7))./( 0.900000.*STATES(:,19)+ 0.100000.*STATES(:,28)));
    ALGEBRAIC(:,58) =  CONSTANTS(:,56).*STATES(:,20).*(STATES(:,1) - ALGEBRAIC(:,57));
    ALGEBRAIC(:,59) =  CONSTANTS(:,57).*(STATES(:,1) - ALGEBRAIC(:,57));
    ALGEBRAIC(:,68) = 1.00000./(1.00000+ 0.124500.*exp((  - 0.100000.*STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11)))+ 0.0365000.*CONSTANTS(:,72).*exp((  - STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))));
    ALGEBRAIC(:,69) = ( (( CONSTANTS(:,66).*ALGEBRAIC(:,68).*1.00000)./(1.00000+power(CONSTANTS(:,67)./STATES(:,19), 1.50000))).*CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,68));
    ALGEBRAIC(:,60) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,28));
    ALGEBRAIC(:,61) =  CONSTANTS(:,58).*power(STATES(:,29), 3.00000).*STATES(:,30).*(STATES(:,1) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,62) =  CONSTANTS(:,59).*STATES(:,31).*STATES(:,32).*(STATES(:,1) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,63) = ( (( 0.293800.*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.000)).*(STATES(:,1) - ALGEBRAIC(:,60)))./(1.00000+exp( 0.0896000.*(STATES(:,1) - ALGEBRAIC(:,60))));
    ALGEBRAIC(:,64) =  CONSTANTS(:,60).*power(STATES(:,33), 2.00000).*(STATES(:,1) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,65) =  CONSTANTS(:,61).*STATES(:,34).*STATES(:,35).*(STATES(:,1) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,66) =  CONSTANTS(:,62).*STATES(:,36).*STATES(:,37).*(STATES(:,1) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,67) =  CONSTANTS(:,63).*STATES(:,38).*(STATES(:,1) -  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.980000.*CONSTANTS(:,7)+ 0.0200000.*CONSTANTS(:,8))./( 0.980000.*STATES(:,28)+ 0.0200000.*STATES(:,19))));
    ALGEBRAIC(:,1) = piecewise({t>=CONSTANTS(:,13)&t<=CONSTANTS(:,14)&(t - CONSTANTS(:,13)) -  floor((t - CONSTANTS(:,13))./CONSTANTS(:,15)).*CONSTANTS(:,15)<=CONSTANTS(:,16), CONSTANTS(:,17) }, 0.00000);
    ALGEBRAIC(:,70) = 0.200000./(1.00000+exp( - (STATES(:,1) - 46.7000)./7.80000));
    ALGEBRAIC(:,71) =  (( CONSTANTS(:,69).*ALGEBRAIC(:,70).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,71))).*(STATES(:,1) - CONSTANTS(:,70));
    ALGEBRAIC(:,72) = arrayfun(@volt_clamp, t);
end

% Compute result of a piecewise function (for stimulation current)
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

function volt = volt_clamp(t)
    if t < 10
        volt = -140;
    elseif (t >=10) && (t <= 80) 
        volt = -60;
    else
        volt = -20;
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_t, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_t = ''; LEGEND_CONSTANTS = '';
    LEGEND_t = strpad('time in component environment (millisecond)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('Cm in component membrane (microF_per_cm2)');
    LEGEND_CONSTANTS(:,2) = strpad('Vmyo in component membrane (microlitre)');
    LEGEND_CONSTANTS(:,3) = strpad('VJSR in component membrane (microlitre)');
    LEGEND_CONSTANTS(:,4) = strpad('VNSR in component membrane (microlitre)');
    LEGEND_CONSTANTS(:,5) = strpad('Vss in component membrane (microlitre)');
    LEGEND_CONSTANTS(:,6) = strpad('Acap in component membrane (cm2)');
    LEGEND_CONSTANTS(:,7) = strpad('Ko in component membrane (micromolar)');
    LEGEND_CONSTANTS(:,8) = strpad('Nao in component membrane (micromolar)');
    LEGEND_CONSTANTS(:,9) = strpad('Cao in component membrane (micromolar)');
    LEGEND_CONSTANTS(:,10) = strpad('R in component membrane (joule_per_mole_kelvin)');
    LEGEND_CONSTANTS(:,11) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,12) = strpad('F in component membrane (coulomb_per_millimole)');
    LEGEND_ALGEBRAIC(:,1) = strpad('i_stim in component membrane (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,47) = strpad('i_CaL in component L_type_calcium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,49) = strpad('i_pCa in component calcium_pump_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,51) = strpad('i_NaCa in component sodium_calcium_exchange_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,55) = strpad('i_Cab in component calcium_background_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,58) = strpad('i_Na in component fast_sodium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,59) = strpad('i_Nab in component sodium_background_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,69) = strpad('i_NaK in component sodium_potassium_pump_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,61) = strpad('i_Kto_f in component fast_transient_outward_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,62) = strpad('i_Kto_s in component slow_transient_outward_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,63) = strpad('i_K1 in component time_independent_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,64) = strpad('i_Ks in component slow_delayed_rectifier_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,65) = strpad('i_Kur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,66) = strpad('i_Kss in component non_inactivating_steady_state_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,71) = strpad('i_ClCa in component calcium_activated_chloride_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,67) = strpad('i_Kr in component rapid_delayed_rectifier_potassium_current (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,13) = strpad('stim_start in component membrane (millisecond)');
    LEGEND_CONSTANTS(:,14) = strpad('stim_end in component membrane (millisecond)');
    LEGEND_CONSTANTS(:,15) = strpad('stim_period in component membrane (millisecond)');
    LEGEND_CONSTANTS(:,16) = strpad('stim_duration in component membrane (millisecond)');
    LEGEND_CONSTANTS(:,17) = strpad('stim_amplitude in component membrane (picoA_per_picoF)');
    LEGEND_STATES(:,2) = strpad('Cai in component calcium_concentration (micromolar)');
    LEGEND_STATES(:,3) = strpad('Cass in component calcium_concentration (micromolar)');
    LEGEND_STATES(:,4) = strpad('CaJSR in component calcium_concentration (micromolar)');
    LEGEND_STATES(:,5) = strpad('CaNSR in component calcium_concentration (micromolar)');
    LEGEND_ALGEBRAIC(:,12) = strpad('Bi in component calcium_concentration (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('Bss in component calcium_concentration (dimensionless)');
    LEGEND_ALGEBRAIC(:,30) = strpad('BJSR in component calcium_concentration (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('CMDN_tot in component calcium_concentration (micromolar)');
    LEGEND_CONSTANTS(:,19) = strpad('CSQN_tot in component calcium_concentration (micromolar)');
    LEGEND_CONSTANTS(:,20) = strpad('Km_CMDN in component calcium_concentration (micromolar)');
    LEGEND_CONSTANTS(:,21) = strpad('Km_CSQN in component calcium_concentration (micromolar)');
    LEGEND_ALGEBRAIC(:,41) = strpad('J_leak in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,34) = strpad('J_rel in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,43) = strpad('J_up in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,37) = strpad('J_tr in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,45) = strpad('J_trpn in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,39) = strpad('J_xfer in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_CONSTANTS(:,22) = strpad('k_plus_htrpn in component calcium_fluxes (per_micromolar_millisecond)');
    LEGEND_CONSTANTS(:,23) = strpad('k_minus_htrpn in component calcium_fluxes (per_millisecond)');
    LEGEND_CONSTANTS(:,24) = strpad('k_plus_ltrpn in component calcium_fluxes (per_micromolar_millisecond)');
    LEGEND_CONSTANTS(:,25) = strpad('k_minus_ltrpn in component calcium_fluxes (per_millisecond)');
    LEGEND_STATES(:,6) = strpad('P_RyR in component calcium_fluxes (dimensionless)');
    LEGEND_CONSTANTS(:,26) = strpad('v1 in component calcium_fluxes (per_millisecond)');
    LEGEND_CONSTANTS(:,27) = strpad('tau_tr in component calcium_fluxes (millisecond)');
    LEGEND_CONSTANTS(:,28) = strpad('v2 in component calcium_fluxes (per_millisecond)');
    LEGEND_CONSTANTS(:,29) = strpad('tau_xfer in component calcium_fluxes (millisecond)');
    LEGEND_CONSTANTS(:,30) = strpad('v3 in component calcium_fluxes (micromolar_per_millisecond)');
    LEGEND_CONSTANTS(:,31) = strpad('Km_up in component calcium_fluxes (micromolar)');
    LEGEND_CONSTANTS(:,32) = strpad('LTRPN_tot in component calcium_buffering (micromolar)');
    LEGEND_CONSTANTS(:,33) = strpad('HTRPN_tot in component calcium_buffering (micromolar)');
    LEGEND_STATES(:,7) = strpad('LTRPN_Ca in component calcium_buffering (micromolar)');
    LEGEND_STATES(:,8) = strpad('HTRPN_Ca in component calcium_buffering (micromolar)');
    LEGEND_CONSTANTS(:,34) = strpad('i_CaL_max in component L_type_calcium_current (picoA_per_picoF)');
    LEGEND_STATES(:,9) = strpad('P_O1 in component ryanodine_receptors (dimensionless)');
    LEGEND_STATES(:,10) = strpad('P_O2 in component ryanodine_receptors (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('P_C1 in component ryanodine_receptors (dimensionless)');
    LEGEND_STATES(:,11) = strpad('P_C2 in component ryanodine_receptors (dimensionless)');
    LEGEND_CONSTANTS(:,35) = strpad('k_plus_a in component ryanodine_receptors (micromolar4_per_millisecond)');
    LEGEND_CONSTANTS(:,36) = strpad('k_minus_a in component ryanodine_receptors (per_millisecond)');
    LEGEND_CONSTANTS(:,37) = strpad('k_plus_b in component ryanodine_receptors (micromolar3_per_millisecond)');
    LEGEND_CONSTANTS(:,38) = strpad('k_minus_b in component ryanodine_receptors (per_millisecond)');
    LEGEND_CONSTANTS(:,39) = strpad('k_plus_c in component ryanodine_receptors (per_millisecond)');
    LEGEND_CONSTANTS(:,40) = strpad('k_minus_c in component ryanodine_receptors (per_millisecond)');
    LEGEND_CONSTANTS(:,41) = strpad('m in component ryanodine_receptors (dimensionless)');
    LEGEND_CONSTANTS(:,42) = strpad('n in component ryanodine_receptors (dimensionless)');
    LEGEND_CONSTANTS(:,43) = strpad('E_CaL in component L_type_calcium_current (millivolt)');
    LEGEND_CONSTANTS(:,44) = strpad('g_CaL in component L_type_calcium_current (milliS_per_microF)');
    LEGEND_STATES(:,12) = strpad('O in component L_type_calcium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,3) = strpad('C1 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,13) = strpad('C2 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,14) = strpad('C3 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,15) = strpad('C4 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,16) = strpad('I1 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,17) = strpad('I2 in component L_type_calcium_current (dimensionless)');
    LEGEND_STATES(:,18) = strpad('I3 in component L_type_calcium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('alpha in component L_type_calcium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,26) = strpad('beta in component L_type_calcium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,31) = strpad('gamma in component L_type_calcium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,35) = strpad('Kpcf in component L_type_calcium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,45) = strpad('Kpcb in component L_type_calcium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,46) = strpad('Kpc_max in component L_type_calcium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,47) = strpad('Kpc_half in component L_type_calcium_current (micromolar)');
    LEGEND_CONSTANTS(:,48) = strpad('i_pCa_max in component calcium_pump_current (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,49) = strpad('Km_pCa in component calcium_pump_current (micromolar)');
    LEGEND_CONSTANTS(:,50) = strpad('k_NaCa in component sodium_calcium_exchange_current (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,51) = strpad('K_mNa in component sodium_calcium_exchange_current (micromolar)');
    LEGEND_CONSTANTS(:,52) = strpad('K_mCa in component sodium_calcium_exchange_current (micromolar)');
    LEGEND_CONSTANTS(:,53) = strpad('k_sat in component sodium_calcium_exchange_current (dimensionless)');
    LEGEND_CONSTANTS(:,54) = strpad('eta in component sodium_calcium_exchange_current (dimensionless)');
    LEGEND_STATES(:,19) = strpad('Nai in component sodium_concentration (micromolar)');
    LEGEND_CONSTANTS(:,55) = strpad('g_Cab in component calcium_background_current (milliS_per_microF)');
    LEGEND_ALGEBRAIC(:,53) = strpad('E_CaN in component calcium_background_current (millivolt)');
    LEGEND_ALGEBRAIC(:,57) = strpad('E_Na in component fast_sodium_current (millivolt)');
    LEGEND_CONSTANTS(:,56) = strpad('g_Na in component fast_sodium_current (milliS_per_microF)');
    LEGEND_STATES(:,20) = strpad('O_Na in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,21) = strpad('C_Na1 in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,22) = strpad('C_Na2 in component fast_sodium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('C_Na3 in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,23) = strpad('I1_Na in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,24) = strpad('I2_Na in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,25) = strpad('IF_Na in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,26) = strpad('IC_Na2 in component fast_sodium_current (dimensionless)');
    LEGEND_STATES(:,27) = strpad('IC_Na3 in component fast_sodium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('alpha_Na11 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,36) = strpad('beta_Na11 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,27) = strpad('alpha_Na12 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,38) = strpad('beta_Na12 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,32) = strpad('alpha_Na13 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,40) = strpad('beta_Na13 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,42) = strpad('alpha_Na3 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,44) = strpad('beta_Na3 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,46) = strpad('alpha_Na2 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,48) = strpad('beta_Na2 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,50) = strpad('alpha_Na4 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,52) = strpad('beta_Na4 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,54) = strpad('alpha_Na5 in component fast_sodium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,56) = strpad('beta_Na5 in component fast_sodium_current (per_millisecond)');
    LEGEND_STATES(:,28) = strpad('Ki in component potassium_concentration (micromolar)');
    LEGEND_CONSTANTS(:,57) = strpad('g_Nab in component sodium_background_current (milliS_per_microF)');
    LEGEND_ALGEBRAIC(:,60) = strpad('E_K in component fast_transient_outward_potassium_current (millivolt)');
    LEGEND_CONSTANTS(:,58) = strpad('g_Kto_f in component fast_transient_outward_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,29) = strpad('ato_f in component fast_transient_outward_potassium_current (dimensionless)');
    LEGEND_STATES(:,30) = strpad('ito_f in component fast_transient_outward_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('alpha_a in component fast_transient_outward_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,15) = strpad('beta_a in component fast_transient_outward_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,6) = strpad('alpha_i in component fast_transient_outward_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,16) = strpad('beta_i in component fast_transient_outward_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,7) = strpad('ass in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('iss in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_CONSTANTS(:,59) = strpad('g_Kto_s in component slow_transient_outward_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,31) = strpad('ato_s in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_STATES(:,32) = strpad('ito_s in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('tau_ta_s in component slow_transient_outward_potassium_current (millisecond)');
    LEGEND_ALGEBRAIC(:,18) = strpad('tau_ti_s in component slow_transient_outward_potassium_current (millisecond)');
    LEGEND_CONSTANTS(:,60) = strpad('g_Ks in component slow_delayed_rectifier_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,33) = strpad('nKs in component slow_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('alpha_n in component slow_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,19) = strpad('beta_n in component slow_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,61) = strpad('g_Kur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,34) = strpad('aur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_STATES(:,35) = strpad('iur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('tau_aur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (millisecond)');
    LEGEND_ALGEBRAIC(:,21) = strpad('tau_iur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (millisecond)');
    LEGEND_CONSTANTS(:,62) = strpad('g_Kss in component non_inactivating_steady_state_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,36) = strpad('aKss in component non_inactivating_steady_state_potassium_current (dimensionless)');
    LEGEND_STATES(:,37) = strpad('iKss in component non_inactivating_steady_state_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('tau_Kss in component non_inactivating_steady_state_potassium_current (millisecond)');
    LEGEND_CONSTANTS(:,63) = strpad('g_Kr in component rapid_delayed_rectifier_potassium_current (milliS_per_microF)');
    LEGEND_STATES(:,38) = strpad('O_K in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_STATES(:,39) = strpad('C_K1 in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_STATES(:,40) = strpad('C_K2 in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,10) = strpad('C_K0 in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_STATES(:,41) = strpad('I_K in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('alpha_a0 in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,28) = strpad('beta_a0 in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,64) = strpad('kb in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,65) = strpad('kf in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,11) = strpad('alpha_a1 in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,24) = strpad('beta_a1 in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,29) = strpad('alpha_i in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_ALGEBRAIC(:,33) = strpad('beta_i in component rapid_delayed_rectifier_potassium_current (per_millisecond)');
    LEGEND_CONSTANTS(:,66) = strpad('i_NaK_max in component sodium_potassium_pump_current (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,67) = strpad('Km_Nai in component sodium_potassium_pump_current (micromolar)');
    LEGEND_CONSTANTS(:,68) = strpad('Km_Ko in component sodium_potassium_pump_current (micromolar)');
    LEGEND_ALGEBRAIC(:,68) = strpad('f_NaK in component sodium_potassium_pump_current (dimensionless)');
    LEGEND_CONSTANTS(:,72) = strpad('sigma in component sodium_potassium_pump_current (dimensionless)');
    LEGEND_CONSTANTS(:,69) = strpad('g_ClCa in component calcium_activated_chloride_current (milliS_per_microF)');
    LEGEND_ALGEBRAIC(:,70) = strpad('O_ClCa in component calcium_activated_chloride_current (dimensionless)');
    LEGEND_CONSTANTS(:,70) = strpad('E_Cl in component calcium_activated_chloride_current (millivolt)');
    LEGEND_CONSTANTS(:,71) = strpad('Km_Cl in component calcium_activated_chloride_current (micromolar)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,2) = strpad('d/dt Cai in component calcium_concentration (micromolar)');
    LEGEND_RATES(:,3) = strpad('d/dt Cass in component calcium_concentration (micromolar)');
    LEGEND_RATES(:,4) = strpad('d/dt CaJSR in component calcium_concentration (micromolar)');
    LEGEND_RATES(:,5) = strpad('d/dt CaNSR in component calcium_concentration (micromolar)');
    LEGEND_RATES(:,6) = strpad('d/dt P_RyR in component calcium_fluxes (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt LTRPN_Ca in component calcium_buffering (micromolar)');
    LEGEND_RATES(:,8) = strpad('d/dt HTRPN_Ca in component calcium_buffering (micromolar)');
    LEGEND_RATES(:,9) = strpad('d/dt P_O1 in component ryanodine_receptors (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt P_O2 in component ryanodine_receptors (dimensionless)');
    LEGEND_RATES(:,11) = strpad('d/dt P_C2 in component ryanodine_receptors (dimensionless)');
    LEGEND_RATES(:,12) = strpad('d/dt O in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,13) = strpad('d/dt C2 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,14) = strpad('d/dt C3 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,15) = strpad('d/dt C4 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,16) = strpad('d/dt I1 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,17) = strpad('d/dt I2 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,18) = strpad('d/dt I3 in component L_type_calcium_current (dimensionless)');
    LEGEND_RATES(:,19) = strpad('d/dt Nai in component sodium_concentration (micromolar)');
    LEGEND_RATES(:,22) = strpad('d/dt C_Na2 in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,21) = strpad('d/dt C_Na1 in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,20) = strpad('d/dt O_Na in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,25) = strpad('d/dt IF_Na in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,23) = strpad('d/dt I1_Na in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,24) = strpad('d/dt I2_Na in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,26) = strpad('d/dt IC_Na2 in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,27) = strpad('d/dt IC_Na3 in component fast_sodium_current (dimensionless)');
    LEGEND_RATES(:,28) = strpad('d/dt Ki in component potassium_concentration (micromolar)');
    LEGEND_RATES(:,29) = strpad('d/dt ato_f in component fast_transient_outward_potassium_current (dimensionless)');
    LEGEND_RATES(:,30) = strpad('d/dt ito_f in component fast_transient_outward_potassium_current (dimensionless)');
    LEGEND_RATES(:,31) = strpad('d/dt ato_s in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_RATES(:,32) = strpad('d/dt ito_s in component slow_transient_outward_potassium_current (dimensionless)');
    LEGEND_RATES(:,33) = strpad('d/dt nKs in component slow_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,34) = strpad('d/dt aur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,35) = strpad('d/dt iur in component ultra_rapidly_activating_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,36) = strpad('d/dt aKss in component non_inactivating_steady_state_potassium_current (dimensionless)');
    LEGEND_RATES(:,37) = strpad('d/dt iKss in component non_inactivating_steady_state_potassium_current (dimensionless)');
    LEGEND_RATES(:,40) = strpad('d/dt C_K2 in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,39) = strpad('d/dt C_K1 in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,38) = strpad('d/dt O_K in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_RATES(:,41) = strpad('d/dt I_K in component rapid_delayed_rectifier_potassium_current (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end
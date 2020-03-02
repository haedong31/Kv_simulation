function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction20101130_yang()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()

    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 1000];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-08, 'AbsTol', 1e-08, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATESODE] = ode15s(@(VOI, STATESODE)computeRates(VOI, STATESODE, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATESODE, ALGEBRAIC] = computeRates(VOI, STATESODE, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATESODE, RATESODE', VOI);

    % Plot state variables against variable of integration
   % [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    STATES = statesconv(STATESODE);
    
    figure(1);
    subplot(4,3,1); plot(VOI, STATES.V); axis tight; xlim([0 300]); ylabel('AP');% V
    subplot(4,3,2); plot(VOI, ALGEBRAIC.ina);  axis tight; xlim([0 300]); ylabel('Ina'); % Ina
    subplot(4,3,3); plot(VOI, ALGEBRAIC.ica);   axis tight; xlim([0 300]); ylabel('Ica');% Ica
    subplot(4,3,4); plot(VOI, ALGEBRAIC.ito);  axis tight; xlim([0 300]); ylabel('Ito');% Ito
    subplot(4,3,5); plot(VOI, ALGEBRAIC.ikr); axis tight; xlim([0 300]); ylabel('Ikr');% Ikr
    subplot(4,3,6); plot(VOI, ALGEBRAIC.inaca); axis tight; xlim([0 300]); ylabel('Inaca');% INaCa
    subplot(4,3,7); plot(VOI, ALGEBRAIC.ikurd); axis tight; xlim([0 300]); ylabel('Ikurd');% Ikurd
    subplot(4,3,8); plot(VOI, ALGEBRAIC.iks); axis tight; xlim([0 300]); ylabel('Iks');% Iks
    subplot(4,3,9); plot(VOI, ALGEBRAIC.inak); axis tight; xlim([0 300]); ylabel('Inak');% INak
    subplot(4,3,10); plot(VOI, ALGEBRAIC.iclca); axis tight; xlim([0 300]);  ylabel('Iclca');% IClCa
    subplot(4,3,11); plot(VOI, ALGEBRAIC.ik1); axis tight; xlim([0 300]); ylabel('Ik1');% Ik1
    subplot(4,3,12); plot(VOI, STATES.cai.*1000); axis tight; xlim([0 300]); ylabel('Cai');% Cai
    
    figure(2);
    ma = ALGEBRAIC.alpha_m./(ALGEBRAIC.alpha_m+ALGEBRAIC.beta_m);
    hi = ALGEBRAIC.alpha_h./(ALGEBRAIC.alpha_h+ALGEBRAIC.beta_h);
    ji =  ALGEBRAIC.alpha_j./(ALGEBRAIC.alpha_j+ALGEBRAIC.beta_j);
    plot(STATES.V,ma.^3,'r'); hold on;
    plot(STATES.V,hi.*ji,'b'); hold off;
    title('Ina activation and inactivation curve');
    
    figure(3)
    plot(STATES.V, ALGEBRAIC.xs_inf.^2,'r');
    title('Iks activation curve');
    
    figure(4)
    plot(STATES.V, ALGEBRAIC.oa_inf.^3,'r'); hold on;
    plot(STATES.V, ALGEBRAIC.oi_inf,'b'); hold off;
    title('Ito activation and inactivation curve');

    
    %     xlabel(LEGEND_VOI);
    %     l = legend(LEGEND_STATES);
    %     set(l,'Interpreter','none');
end


function [INIT_STATES, CONSTANTS] = initConsts()
    CONSTANTS = []; 
    STATES = [];     
    
    CONSTANTS.R = 8.3143;
    CONSTANTS.T = 310.0;
    CONSTANTS.F = 96.4867;
    CONSTANTS.Cm = 1.0; %Cm=100
    CONSTANTS.vi = 13668.0;%vi = 13668.0
    CONSTANTS.v_up = 1109.52;
    CONSTANTS.v_rel = 96.48;
    CONSTANTS.ko = 5.4;
    CONSTANTS.nao = 140.0;
    CONSTANTS.cao = 2;%cao = 1.8
    CONSTANTS.clo = 132.0;
    CONSTANTS.gna = 8.8; %gna = 7.8
    CONSTANTS.gk1 = 0.15;
    CONSTANTS.gto = 0.19824;
    CONSTANTS.gkr = 0.06984; %gkr = 0.06984
    CONSTANTS.gks = 0.0661; %gks = 0.0561
    CONSTANTS.gca = 0.09; %gca = 0.24
    CONSTANTS.gclca = 30000000000; %gclca = 0.3
    CONSTANTS.g_b_ca = 0.00113; %g_b_ca = 0.00113
    CONSTANTS.g_b_na = 0.000674;
    CONSTANTS.inakmax = 0.6; %inakmax = 0.6;
    CONSTANTS.inacamax = 1600.0;%inacamax = 1600.0;
    CONSTANTS.i_p_ca_max = 0.275; %i_p_ca_max = 0.275
    CONSTANTS.i_up_max = 0.005; %i_up_max = 0.005   
    CONSTANTS.stim_start = 10; %stim_start in component membrane (millisecond)
    CONSTANTS.stim_period = 1000; %stim_period in component membrane (millisecond)
    CONSTANTS.stim_dur = 1; %stim_duration in component membrane (millisecond)
    CONSTANTS.istim = -80.0; %stim_amplitude in component membrane (picoA_per_picoF) 2000
    
    CONSTANTS.tau_fca = 2.0;
    CONSTANTS.tau_qca = 2.0;
    CONSTANTS.tau_tr = 180.0;
    CONSTANTS.tau_u = 8.00000;    
    CONSTANTS.kmnai = 10.0;
    CONSTANTS.kmko = 1.5;    
    CONSTANTS.kmna = 87.5;%kmna = 87.5
    CONSTANTS.kmca = 1.38;
    CONSTANTS.ksat = 0.1; %ksat = 0.1;
    CONSTANTS.k_rel = 30.0;
    CONSTANTS.k_up = 0.00092;  %k_up = 0.00092  
    CONSTANTS.ca_up_max = 15.0; %ca_up_max = 15.0;
    CONSTANTS.cmdn_max = 0.045;%cmdn_max = 0.045
    CONSTANTS.trpn_max = 0.35;%trpn_max = 0.35
    CONSTANTS.csqn_max = 10.0;    
    CONSTANTS.sigma =  (1.00000./7.00000).*((exp((CONSTANTS.nao./67.3000))) - 1.00000);
    
    STATES.V = -83.53;
    STATES.m = 0.001972;
    STATES.h = 0.9791;
    STATES.j = 0.9869;
    STATES.d = 0.000004757;
    STATES.f = 0.9999;
    STATES.xr = 0.0000007433;
    STATES.xs = 0.01791;
    STATES.nai = 11.75;
    STATES.cai = 0.0001024;%cai = 0.0001024
    STATES.ki = 138.4;
    STATES.cli = 29.26;
    STATES.ca_up = 1.502;
    STATES.ca_rel = 1.502;
    STATES.oa = 0.07164;
    STATES.oi = 0.9980;
    STATES.qca=0;
    STATES.ua = 0.05869;
    STATES.ui = 0.9987;
    STATES.ca_cmdn = 0.001856;
    STATES.ca_trpn = 0.007022;
    STATES.ca_csqn = 6.432;
    STATES.fca = 0.7484;
    STATES.u = 0.0;
    STATES.v = 1.0;
    STATES.w = 0.9993;
    
    INIT_STATES = statesconv(STATES);
end

function [RATESODE, ALGEBRAIC] = computeRates(VOI, STATESODE, CONSTANTS)

    statesSize = size(STATESODE);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATESODE = STATESODE';
    else
        statesRowCount = statesSize(1);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    
    STATES = statesconv(STATESODE);

    ALGEBRAIC.alpha_m =  0.320000.*((STATES.V+47.1300)./(1.00000 - (exp(( -0.100000.*(STATES.V+47.1300))))));%A12
    ALGEBRAIC.beta_m =  0.0800000.*(exp((STATES.V./-11.0000))); %A13
    RATES.dmdt =  ALGEBRAIC.alpha_m.*(1.00000 - STATES.m) -  ALGEBRAIC.beta_m.*STATES.m;
    
    ALGEBRAIC.alpha_h = piecewise({STATES.V<-40.0000,  0.135000.*(exp(((STATES.V+80.0000)./-6.80000))) }, 0.00000);%A14
    ALGEBRAIC.beta_h = piecewise({STATES.V<-40.0000,  3.56000.*(exp(( 0.0790000.*STATES.V)))+...
        310000..*(exp(( 0.350000.*STATES.V))) }, 1.00000./( 0.130000.*(1.00000+(exp(((STATES.V+10.6600)./-11.1000)))))); %A15
    RATES.dhdt =  ALGEBRAIC.alpha_h.*(1.00000 - STATES.h) -  ALGEBRAIC.beta_h.*STATES.h;
    
    ALGEBRAIC.alpha_j = piecewise({STATES.V<-40.0000,  (( -127140..*(exp(( 0.244400.*STATES.V))) -  ...
        3.47400e-05.*(exp(( -0.0439100.*STATES.V))))./(1.00000+(exp(( 0.311000.*(STATES.V+79.2300)))))).*(STATES.V+37.7800) }, 0.00000); %A16
    ALGEBRAIC.beta_j = piecewise({STATES.V<-40.0000, ( 0.121200.*(exp(( -0.0105200.*STATES.V))))./(1.00000+...
        (exp(( -0.137800.*(STATES.V+40.1400))))) }, ( 0.300000.*(exp(( -2.53500e-07.*STATES.V))))./(1.00000+(exp(( -0.100000.*(STATES.V+32.0000)))))); %A17
    RATES.djdt =  ALGEBRAIC.alpha_j.*(1.00000 - STATES.j) -  ALGEBRAIC.beta_j.*STATES.j;
    
    ALGEBRAIC.tau_d = (1.00000 - (exp(((STATES.V+10.0000)./-6.24000))))./( 0.0350000.*(STATES.V+...
        10.0000).*(1.00000+(exp(((STATES.V+10.0000)./-6.24000))))); %A50
    ALGEBRAIC.dinf = (1.00000+(exp(((STATES.V+10.0000)./-6.00000)))) .^ -1.00000; %A51
    RATES.dddt = (ALGEBRAIC.dinf - STATES.d)./ALGEBRAIC.tau_d;    
    
    ALGEBRAIC.tau_f =  400.000.*((1.00000+ 4.50000.*(exp(( -0.000700000.*((STATES.V - 9.00000) .^ 2.00000))))) .^ -1.00000);%A52
    ALGEBRAIC.f_inf = (1.00000+(exp(((STATES.V+24.6000)./6.20000)))) .^ -1.00000;%A53
    RATES.dfdt = (ALGEBRAIC.f_inf - STATES.f)./ALGEBRAIC.tau_f;
    ALGEBRAIC.fca_inf = 0.290000+ 0.800000.*((1.00000+(exp(((STATES.cai - 0.000120000)./6.00000e-05)))) .^ -1.00000); %A54
    RATES.dfcadt = (ALGEBRAIC.fca_inf - STATES.fca)./CONSTANTS.tau_fca;
        
    ALGEBRAIC.alpha_oa =  0.650000.*(((exp(((STATES.V+18.0000)./-8.50000)))+(exp(((STATES.V - 16.0000)./-59.0000)))) .^ -1.00000);%A21
    ALGEBRAIC.beta_oa =  1.20000.*((2.20000+(exp(((STATES.V+75.0000)./18.0000)))) .^ -1.00000);%A22
    ALGEBRAIC.tau_oa = (ALGEBRAIC.alpha_oa+ALGEBRAIC.beta_oa) .^ -1.00000;%A23
    ALGEBRAIC.oa_inf = (1.00000+(exp(((STATES.V+0.500000)./-10.5000)))) .^ ( - (1.00000./3.00000));%A24    
    RATES.doadt = (ALGEBRAIC.oa_inf - STATES.oa)./ALGEBRAIC.tau_oa;
    
    ALGEBRAIC.alpha_oi = (6.20000+(exp(((STATES.V+105.200)./9.85000)))) .^ -1.00000;%A25
    ALGEBRAIC.beta_oi = (7.54000+(exp(((STATES.V - 8.89000)./-12.8700)))) .^ -1.00000;%A26
    ALGEBRAIC.tau_oi = (ALGEBRAIC.alpha_oi+ALGEBRAIC.beta_oi) .^ -1.00000;%A27
    ALGEBRAIC.oi_inf = (1.00000+(exp(((STATES.V+43.3770)./6.45000)))) .^ -1.00000;%A28
    RATES.doidt = (ALGEBRAIC.oi_inf - STATES.oi)./ALGEBRAIC.tau_oi;    
    
    ALGEBRAIC.alpha_ua =  1.47000.*(((exp(((STATES.V+33.2000)./-30.6300)))+(exp(((STATES.V - 27.6000)./-30.6500)))) .^ -1.00000);%A31
    ALGEBRAIC.beta_ua =  0.420000.*(((exp(((STATES.V+26.6400)./2.49000)))+(exp(((STATES.V+44.4100)./20.3600)))) .^ -1.00000);%A32
    ALGEBRAIC.tau_ua = (ALGEBRAIC.alpha_ua+ALGEBRAIC.beta_ua) .^ -1.00000;%A33
    ALGEBRAIC.ua_inf = (1.00000+(exp(((STATES.V+2.81000)./-9.49000)))) .^ ( - (1.00000./3.00000)); %A34
    RATES.duadt = (ALGEBRAIC.ua_inf - STATES.ua)./ALGEBRAIC.tau_ua;
    
    ALGEBRAIC.alpha_ui = (21.0000+(exp(((STATES.V - 185.000)./-28.0000)))) .^ -1.00000;%A35
    ALGEBRAIC.beta_ui = exp(((STATES.V - 158.000)./16.0000));%A36
    ALGEBRAIC.tau_ui = (ALGEBRAIC.alpha_ui+ALGEBRAIC.beta_ui) .^ -1.00000;%A37
    ALGEBRAIC.ui_inf = (1.00000+(exp(((STATES.V - 99.4500)./27.4800)))) .^ -1.00000;%A38
    RATES.duidt = (ALGEBRAIC.ui_inf - STATES.ui)./ALGEBRAIC.tau_ui;    
    
    ALGEBRAIC.alpha_xr =  0.0400000.*((STATES.V - 248.000)./(1.00000 - (exp(((STATES.V - 248.000)./-28.0000))))); %A40
    ALGEBRAIC.beta_xr =  0.0280000.*((STATES.V+163.000)./((exp(((STATES.V+163.000)./21.0000))) - 1.00000));%A41
    ALGEBRAIC.tau_xr = (ALGEBRAIC.alpha_xr+ALGEBRAIC.beta_xr) .^ -1.00000;%A42
    ALGEBRAIC.xr_inf = (1.00000+(exp(((STATES.V+7.65400)./-5.37700)))) .^ -1.00000;%A43
    RATES.dxrdt = (ALGEBRAIC.xr_inf - STATES.xr)./ALGEBRAIC.tau_xr;
    
    ALGEBRAIC.alpha_xs =  1.00000e-05.*((STATES.V+28.5000)./(1.00000 - (exp(((STATES.V+28.5000)./-115.000)))));%A45
    ALGEBRAIC.beta_xs =  0.000230000.*((STATES.V+28.5000)./((exp(((STATES.V+28.5000)./3.30000))) - 1.00000));%A46
    ALGEBRAIC.tau_xs = (ALGEBRAIC.alpha_xs+ALGEBRAIC.beta_xs) .^ -1.00000;%A47
    ALGEBRAIC.xs_inf = (1.00000+(exp(((STATES.V - 13.0000)./-12.0000)))) .^ -0.500000;%A48
    RATES.dxsdt = (ALGEBRAIC.xs_inf - STATES.xs)./ALGEBRAIC.tau_xs;    
    
    ALGEBRAIC.Ek =  (( CONSTANTS.R.*CONSTANTS.T)./CONSTANTS.F).*(log((CONSTANTS.ko./STATES.ki)));%A10
    ALGEBRAIC.ik1 = ( CONSTANTS.gk1.*(STATES.V - ALGEBRAIC.Ek))./(1.00000+(exp(( 0.0700000.*(STATES.V+80.0000))))); %A19
    ALGEBRAIC.ito =  CONSTANTS.gto.*(STATES.oa .^ 3.00000).*STATES.oi.*(STATES.V - ALGEBRAIC.Ek); %A20
    
    ALGEBRAIC.gkurd = 0.00855000+0.0779000./(1.00000+(exp(((STATES.V+11.0000)./-16.0000))));  %A30  
    ALGEBRAIC.ikurd =  ALGEBRAIC.gkurd.*(STATES.ua .^ 3.00000).*STATES.ui.*(STATES.V - ALGEBRAIC.Ek); %A29    
    
    ALGEBRAIC.ikr =  CONSTANTS.gkr.*STATES.xr.*(0.0700000+0.580000./(1.00000+(exp(((STATES.V+...
        15.0000)./22.4000))))).*(STATES.V - ALGEBRAIC.Ek); %A39
    ALGEBRAIC.iks =  CONSTANTS.gks.*(STATES.xs .^ 2.00000).*(STATES.V - ALGEBRAIC.Ek); %A44
    
    ALGEBRAIC.fnak = (1.00000+ 0.124500.*(exp(( -0.100000.*(( CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))))...
        + 0.0365000.*CONSTANTS.sigma.*(exp(( - (( CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))))) .^ -1.00000; %A58
    ALGEBRAIC.inak =  CONSTANTS.inakmax.*ALGEBRAIC.fnak.*(1.00000./(1.00000+...
        ((CONSTANTS.kmnai./STATES.nai) .^ 1.50000))).*(CONSTANTS.ko./(CONSTANTS.ko+CONSTANTS.kmko)); %A57
    
    RATES.dkidt = ( 2.00000.*ALGEBRAIC.inak - (ALGEBRAIC.ik1+ALGEBRAIC.ito+ALGEBRAIC.ikurd+ALGEBRAIC.ikr+...
        ALGEBRAIC.iks))./( CONSTANTS.vi.*CONSTANTS.F); % A4
    ALGEBRAIC.Ena =  (( CONSTANTS.R.*CONSTANTS.T)./CONSTANTS.F).*(log((CONSTANTS.nao./STATES.nai)));%A10
    ALGEBRAIC.ina =  CONSTANTS.gna.*(STATES.m .^ 3.00000).*STATES.h.*STATES.j.*(STATES.V - ALGEBRAIC.Ena); %A11
    
    ALGEBRAIC.inaca = ( CONSTANTS.inacamax.*( (exp((( 0.350000.*CONSTANTS.F.*STATES.V)./...
        ( CONSTANTS.R.*CONSTANTS.T)))).*(STATES.nai .^ 3.00000).*CONSTANTS.cao -  ...
        (exp((( -0.650000.*CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))).*(CONSTANTS.nao .^ 3.00000).*STATES.cai))...
        ./( ((CONSTANTS.kmna .^ 3.00000)+(CONSTANTS.nao .^ 3.00000)).*(CONSTANTS.kmca+CONSTANTS.cao).*(100.00000+ ...
        CONSTANTS.ksat.*(exp((( -0.650000.*STATES.V.*CONSTANTS.F)./( CONSTANTS.R.*CONSTANTS.T))))));%A60
    
    ALGEBRAIC.ibna =  CONSTANTS.g_b_na.*(STATES.V - ALGEBRAIC.Ena); %A62
    RATES.dnaidt = ( -3.00000.*ALGEBRAIC.inak - ( 3.00000.*ALGEBRAIC.inaca+ALGEBRAIC.ibna+...
        ALGEBRAIC.ina))./( CONSTANTS.vi.*CONSTANTS.F); %A3
    
    ALGEBRAIC.ica =  CONSTANTS.gca.*STATES.d.*STATES.f.*STATES.fca.*(STATES.V - 65.0000);%A49
    ALGEBRAIC.ipca =  CONSTANTS.i_p_ca_max.*(STATES.cai./(0.000500000+STATES.cai));%A63
    ALGEBRAIC.Eca =  (( CONSTANTS.R.*CONSTANTS.T)./( 2.00000.*CONSTANTS.F)).*(log((CONSTANTS.cao./STATES.cai)));  %A10
    ALGEBRAIC.ibca =  CONSTANTS.g_b_ca.*(STATES.V - ALGEBRAIC.Eca); %A61
    
    RATES.dca_csqndt =  0.480000.*STATES.ca_rel.*(1.00000 - STATES.ca_csqn./CONSTANTS.csqn_max) -  ...
        0.400000.*(STATES.ca_csqn./CONSTANTS.csqn_max); %A76
    
    ALGEBRAIC.i_rel =  CONSTANTS.k_rel.*(STATES.u .^ 2.00000).*STATES.v.*STATES.w.*(STATES.ca_rel - STATES.cai); %A64
    ALGEBRAIC.Fn =  1.00000e-12.*CONSTANTS.v_rel.*ALGEBRAIC.i_rel -  5.00000e-13.*( (1.00000./( 2.00000.*CONSTANTS.F))...
        .*ALGEBRAIC.ica -  (1.00000./( 5.00000.*CONSTANTS.F)).*ALGEBRAIC.inaca); %A70
    ALGEBRAIC.u_inf = (1.00000+(exp(((ALGEBRAIC.Fn - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000; %A65
    RATES.dudt = (ALGEBRAIC.u_inf - STATES.u)./CONSTANTS.tau_u;
    ALGEBRAIC.tau_v = 1.91000+ 2.09000.*((1.00000+(exp(((ALGEBRAIC.Fn - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000);%A66
    ALGEBRAIC.v_inf = 1.00000 - ((1.00000+(exp(((ALGEBRAIC.Fn - 6.83500e-14)./-1.36700e-15)))) .^ -1.00000);%A67
    RATES.djsrvdt = (ALGEBRAIC.v_inf - STATES.v)./ALGEBRAIC.tau_v;
    ALGEBRAIC.tau_w = (6.00000 -  6.00000.*(exp(((STATES.V - 7.90000)./-5.00000))))./( (1.00000+ ...
        0.300000.*(exp(((STATES.V - 7.90000)./-5.00000)))).*(STATES.V - 7.90000)); %A68
    ALGEBRAIC.w_inf = 1.00000 - ((1.00000+(exp(((STATES.V - 40.0000)./-17.0000)))) .^ -1.00000); %A69
    RATES.dwdt = (ALGEBRAIC.w_inf - STATES.w)./ALGEBRAIC.tau_w;
    
    ALGEBRAIC.Ecl =  (( CONSTANTS.R.*CONSTANTS.T)./( -1.00000.*CONSTANTS.F)).*(log((CONSTANTS.clo./STATES.cli)));%A10
    ALGEBRAIC.iclca =  CONSTANTS.gclca.*STATES.qca.*(STATES.V - ALGEBRAIC.Ecl); %A55
    ALGEBRAIC.qca_inf = 1-(1+(ALGEBRAIC.Fn./1.00000e-10).^3).^(-1); %A56
    RATES.dqcadt = (ALGEBRAIC.qca_inf - STATES.qca)./CONSTANTS.tau_qca;
    RATES.dclidt = ALGEBRAIC.iclca./( CONSTANTS.vi.*CONSTANTS.F); % A5
    
    RATES.dca_cmdndt =  200.000.*STATES.cai.*(1.00000 - STATES.ca_cmdn./CONSTANTS.cmdn_max) -  ...
        0.476000.*(STATES.ca_cmdn./CONSTANTS.cmdn_max); %A74
    ALGEBRAIC.i_tr = (STATES.ca_up - STATES.ca_rel)./CONSTANTS.tau_tr; %A71
    ALGEBRAIC.i_up = CONSTANTS.i_up_max./(1.00000+CONSTANTS.k_up./STATES.cai); %A72
    ALGEBRAIC.i_up_leak =  CONSTANTS.i_up_max.*(STATES.ca_up./CONSTANTS.ca_up_max); %A73
    RATES.dcaupdt = ALGEBRAIC.i_up - (ALGEBRAIC.i_up_leak+ ALGEBRAIC.i_tr.*(CONSTANTS.v_rel./CONSTANTS.v_up)); % A7
    
    ALGEBRAIC.j_ca_csqn = RATES.dca_csqndt;
    RATES.dcareldt = ALGEBRAIC.i_tr - (ALGEBRAIC.i_rel+ 31.0000.*ALGEBRAIC.j_ca_csqn); % A8
    
    RATES.dca_trpndt =  78.4000.*STATES.cai.*(1.00000 - STATES.ca_trpn./CONSTANTS.trpn_max) ...
        -  0.392000.*(STATES.ca_trpn./CONSTANTS.trpn_max);%A75
    ALGEBRAIC.j_ca_cmdn = RATES.dca_cmdndt;
    ALGEBRAIC.j_ca_trpn = RATES.dca_trpndt;
    
    RATES.dcaidt = (( 2.00000.*ALGEBRAIC.inaca - (ALGEBRAIC.ipca+ALGEBRAIC.ica+ALGEBRAIC.ibca))...
        ./( 2.00000.*CONSTANTS.vi.*CONSTANTS.F)+( CONSTANTS.v_up.*(ALGEBRAIC.i_up_leak - ALGEBRAIC.i_up)...
        + ALGEBRAIC.i_rel.*CONSTANTS.v_rel)./CONSTANTS.vi) - ( CONSTANTS.trpn_max.*ALGEBRAIC.j_ca_trpn+ ...
        CONSTANTS.cmdn_max.*ALGEBRAIC.j_ca_cmdn); % A6
    
     STIM = piecewise({VOI -  (floor((VOI./CONSTANTS.stim_period))).*CONSTANTS.stim_period>=CONSTANTS.stim_start& VOI -  ...
        (floor((VOI./CONSTANTS.stim_period))).*CONSTANTS.stim_period<=CONSTANTS.stim_start+CONSTANTS.stim_dur,  CONSTANTS.istim}, 0.00000); % A2
    
    RATES.dvdt =  - (ALGEBRAIC.ina+ALGEBRAIC.ik1+ALGEBRAIC.ito+ALGEBRAIC.ikurd+ALGEBRAIC.ikr+ALGEBRAIC.iks+ALGEBRAIC.ica...
        +ALGEBRAIC.iclca+ALGEBRAIC.ipca+ALGEBRAIC.inaca+ALGEBRAIC.inak+ALGEBRAIC.ibna+ALGEBRAIC.ibca+STIM)./CONSTANTS.Cm; % A1
    
    RATESODE = [RATES.dvdt, RATES.dmdt, RATES.dhdt, RATES.djdt, RATES.dddt, RATES.dfdt, RATES.dxrdt, RATES.dxsdt,  ...
        RATES.dnaidt, RATES.dcaidt, RATES.dkidt, RATES.dclidt, RATES.dcaupdt, RATES.dcareldt, ...
        RATES.doadt, RATES.doidt, RATES.dqcadt, RATES.duadt, RATES.duidt, RATES.dfcadt,  RATES.dca_cmdndt, ...
        RATES.dca_trpndt, RATES.dca_csqndt, RATES.dudt, RATES.djsrvdt, RATES.dwdt];
    RATESODE = RATESODE';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATESODE, RATESODE, VOI)
    STATES = statesconv(STATESODE);
    
    RATES.dca_csqndt = RATESODE(:,23);
    RATES.dca_cmdndt = RATESODE(:,21);
    RATES.dca_trpndt = RATESODE(:,22);    
    
    ALGEBRAIC.fca_inf = 0.290000+ 0.800000.*((1.00000+(exp(((STATES.cai - 0.000120000)./6.00000e-05)))) .^ -1.00000);
    ALGEBRAIC.alpha_m =  0.320000.*((STATES.V+47.1300)./(1.00000 - (exp(( -0.100000.*(STATES.V+47.1300))))));
    ALGEBRAIC.beta_m =  0.0800000.*(exp((STATES.V./-11.0000)));
    ALGEBRAIC.alpha_h = piecewise({STATES.V<-40.0000,  0.135000.*(exp(((STATES.V+80.0000)./-6.80000))) }, 0.00000);
    ALGEBRAIC.beta_h = piecewise({STATES.V<-40.0000,  3.56000.*(exp(( 0.0790000.*STATES.V)))+ 310000..*(exp(( 0.350000.*STATES.V))) }, 1.00000./( 0.130000.*(1.00000+(exp(((STATES.V+10.6600)./-11.1000))))));
    ALGEBRAIC.alpha_j = piecewise({STATES.V<-40.0000,  (( -127140..*(exp(( 0.244400.*STATES.V))) -  3.47400e-05.*(exp(( -0.0439100.*STATES.V))))./(1.00000+(exp(( 0.311000.*(STATES.V+79.2300)))))).*(STATES.V+37.7800) }, 0.00000);
    ALGEBRAIC.beta_j = piecewise({STATES.V<-40.0000, ( 0.121200.*(exp(( -0.0105200.*STATES.V))))./(1.00000+(exp(( -0.137800.*(STATES.V+40.1400))))) }, ( 0.300000.*(exp(( -2.53500e-07.*STATES.V))))./(1.00000+(exp(( -0.100000.*(STATES.V+32.0000))))));
    ALGEBRAIC.dinf = (1.00000+(exp(((STATES.V+10.0000)./-6.00000)))) .^ -1.00000;
    ALGEBRAIC.tau_d = (1.00000 - (exp(((STATES.V+10.0000)./-6.24000))))./( 0.0350000.*(STATES.V+10.0000).*(1.00000+(exp(((STATES.V+10.0000)./-6.24000)))));
    ALGEBRAIC.f_inf = (1.00000+(exp(((STATES.V+24.6000)./6.20000)))) .^ -1.00000;
    ALGEBRAIC.tau_f =  400.000.*((1.00000+ 4.50000.*(exp(( -0.000700000.*((STATES.V - 9.00000) .^ 2.00000))))) .^ -1.00000);
    ALGEBRAIC.tau_w = (6.00000 -  6.00000.*(exp(((STATES.V - 7.90000)./-5.00000))))./( (1.00000+ 0.300000.*(exp(((STATES.V - 7.90000)./-5.00000)))).*(STATES.V - 7.90000));
    ALGEBRAIC.w_inf = 1.00000 - ((1.00000+(exp(((STATES.V - 40.0000)./-17.0000)))) .^ -1.00000);
    ALGEBRAIC.alpha_oa =  0.650000.*(((exp(((STATES.V+18.0000)./-8.50000)))+(exp(((STATES.V - 16.0000)./-59.0000)))) .^ -1.00000);
    ALGEBRAIC.beta_oa =  1.20000.*((2.20000+(exp(((STATES.V+75.0000)./18.0000)))) .^ -1.00000);
    ALGEBRAIC.tau_oa = (ALGEBRAIC.alpha_oa+ALGEBRAIC.beta_oa) .^ -1.00000;
    ALGEBRAIC.oa_inf = (1.00000+(exp(((STATES.V+0.500000)./-10.5000)))) .^ ( - (1.00000./3.00000));
    ALGEBRAIC.alpha_oi = (6.20000+(exp(((STATES.V+105.200)./9.85000)))) .^ -1.00000;
    ALGEBRAIC.beta_oi = (7.54000+(exp(((STATES.V - 8.89000)./-12.8700)))) .^ -1.00000;
    ALGEBRAIC.tau_oi = (ALGEBRAIC.alpha_oi+ALGEBRAIC.beta_oi) .^ -1.00000;
    ALGEBRAIC.oi_inf = (1.00000+(exp(((STATES.V+43.3770)./6.45000)))) .^ -1.00000;
    ALGEBRAIC.alpha_ua =  1.47000.*(((exp(((STATES.V+33.2000)./-30.6300)))+(exp(((STATES.V - 27.6000)./-30.6500)))) .^ -1.00000);
    ALGEBRAIC.beta_ua =  0.420000.*(((exp(((STATES.V+26.6400)./2.49000)))+(exp(((STATES.V+44.4100)./20.3600)))) .^ -1.00000);
    ALGEBRAIC.tau_ua = (ALGEBRAIC.alpha_ua+ALGEBRAIC.beta_ua) .^ -1.00000;
    ALGEBRAIC.ua_inf = (1.00000+(exp(((STATES.V+2.81000)./-9.49000)))) .^ ( - (1.00000./3.00000));
    ALGEBRAIC.alpha_ui = (21.0000+(exp(((STATES.V - 185.000)./-28.0000)))) .^ -1.00000;
    ALGEBRAIC.beta_ui = exp(((STATES.V - 158.000)./16.0000));
    ALGEBRAIC.tau_ui = (ALGEBRAIC.alpha_ui+ALGEBRAIC.beta_ui) .^ -1.00000;
    ALGEBRAIC.ui_inf = (1.00000+(exp(((STATES.V - 99.4500)./27.4800)))) .^ -1.00000;
    ALGEBRAIC.alpha_xr =  0.0400000.*((STATES.V - 248.000)./(1.00000 - (exp(((STATES.V - 248.000)./-28.0000)))));
    ALGEBRAIC.beta_xr =  0.0280000.*((STATES.V+163.000)./((exp(((STATES.V+163.000)./21.0000))) - 1.00000));
    ALGEBRAIC.tau_xr = (ALGEBRAIC.alpha_xr+ALGEBRAIC.beta_xr) .^ -1.00000;
    ALGEBRAIC.xr_inf = (1.00000+(exp(((STATES.V+7.65400)./-5.37700)))) .^ -1.00000;
    ALGEBRAIC.alpha_xs =  1.00000e-05.*((STATES.V+28.5000)./(1.00000 - (exp(((STATES.V+28.5000)./-115.000)))));
    ALGEBRAIC.beta_xs =  0.000230000.*((STATES.V+28.5000)./((exp(((STATES.V+28.5000)./3.30000))) - 1.00000));
    ALGEBRAIC.tau_xs = (ALGEBRAIC.alpha_xs+ALGEBRAIC.beta_xs) .^ -1.00000;
    ALGEBRAIC.xs_inf = (1.00000+(exp(((STATES.V - 13.0000)./-12.0000)))) .^ -0.500000;
    ALGEBRAIC.Ecl =  (( CONSTANTS.R.*CONSTANTS.T)./( -1.00000.*CONSTANTS.F)).*(log((CONSTANTS.clo./STATES.cli)));
    ALGEBRAIC.iclca =  CONSTANTS.gclca.*STATES.qca.*(STATES.V - ALGEBRAIC.Ecl); %A55
    ALGEBRAIC.Ek =  (( CONSTANTS.R.*CONSTANTS.T)./CONSTANTS.F).*(log((CONSTANTS.ko./STATES.ki)));
    ALGEBRAIC.ik1 = ( CONSTANTS.gk1.*(STATES.V - ALGEBRAIC.Ek))./(1.00000+(exp(( 0.0700000.*(STATES.V+80.0000)))));
    ALGEBRAIC.ito =  CONSTANTS.gto.*(STATES.oa .^ 3.00000).*STATES.oi.*(STATES.V - ALGEBRAIC.Ek);
    ALGEBRAIC.gkurd = 0.00855000+0.0779000./(1.00000+(exp(((STATES.V+11.0000)./-16.0000))));
    ALGEBRAIC.ikurd =  ALGEBRAIC.gkurd.*(STATES.ua .^ 3.00000).*STATES.ui.*(STATES.V - ALGEBRAIC.Ek);
    ALGEBRAIC.ikr =  CONSTANTS.gkr.*STATES.xr.*(0.0700000+0.580000./(1.00000+(exp(((STATES.V+15.0000)./22.4000))))).*(STATES.V - ALGEBRAIC.Ek);
    ALGEBRAIC.iks =  CONSTANTS.gks.*(STATES.xs .^ 2.00000).*(STATES.V - ALGEBRAIC.Ek);
    ALGEBRAIC.fnak = (1.00000+ 0.124500.*(exp(( -0.100000.*(( CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))))+ 0.0365000.*CONSTANTS.sigma.*(exp(( - (( CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))))) .^ -1.00000;
    ALGEBRAIC.inak =  CONSTANTS.inakmax.*ALGEBRAIC.fnak.*(1.00000./(1.00000+((CONSTANTS.kmnai./STATES.nai) .^ 1.50000))).*(CONSTANTS.ko./(CONSTANTS.ko+CONSTANTS.kmko));
    ALGEBRAIC.Ena =  (( CONSTANTS.R.*CONSTANTS.T)./CONSTANTS.F).*(log((CONSTANTS.nao./STATES.nai)));
    ALGEBRAIC.ina =  CONSTANTS.gna.*(STATES.m .^ 3.00000).*STATES.h.*STATES.j.*(STATES.V - ALGEBRAIC.Ena);
    ALGEBRAIC.inaca = ( CONSTANTS.inacamax.*( (exp((( 0.350000.*CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))).*(STATES.nai .^ 3.00000).*CONSTANTS.cao -  (exp((( -0.650000.*CONSTANTS.F.*STATES.V)./( CONSTANTS.R.*CONSTANTS.T)))).*(CONSTANTS.nao .^ 3.00000).*STATES.cai))./( ((CONSTANTS.kmna .^ 3.00000)+(CONSTANTS.nao .^ 3.00000)).*(CONSTANTS.kmca+CONSTANTS.cao).*(1.00000+ CONSTANTS.ksat.*(exp((( -0.650000.*STATES.V.*CONSTANTS.F)./( CONSTANTS.R.*CONSTANTS.T))))));
    ALGEBRAIC.ibna =  CONSTANTS.g_b_na.*(STATES.V - ALGEBRAIC.Ena);
    ALGEBRAIC.ica =  CONSTANTS.gca.*STATES.d.*STATES.f.*STATES.fca.*(STATES.V - 65.0000);
    ALGEBRAIC.ipca =  CONSTANTS.i_p_ca_max.*(STATES.cai./(0.000500000+STATES.cai));
    ALGEBRAIC.Eca =  (( CONSTANTS.R.*CONSTANTS.T)./( 2.00000.*CONSTANTS.F)).*(log((CONSTANTS.cao./STATES.cai)));
    ALGEBRAIC.ibca =  CONSTANTS.g_b_ca.*(STATES.V - ALGEBRAIC.Eca);
    ALGEBRAIC.i_rel =  CONSTANTS.k_rel.*(STATES.u .^ 2.00000).*STATES.v.*STATES.w.*(STATES.ca_rel - STATES.cai);
    ALGEBRAIC.Fn =  1.00000e-12.*CONSTANTS.v_rel.*ALGEBRAIC.i_rel -  5.00000e-13.*( (1.00000./( 2.00000.*CONSTANTS.F)).*ALGEBRAIC.ica -  (1.00000./( 5.00000.*CONSTANTS.F)).*ALGEBRAIC.inaca);
    ALGEBRAIC.u_inf = (1.00000+(exp(((ALGEBRAIC.Fn - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000;
    ALGEBRAIC.tau_v = 1.91000+ 2.09000.*((1.00000+(exp(((ALGEBRAIC.Fn - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000);
    ALGEBRAIC.v_inf = 1.00000 - ((1.00000+(exp(((ALGEBRAIC.Fn - 6.83500e-14)./-1.36700e-15)))) .^ -1.00000);
    ALGEBRAIC.i_tr = (STATES.ca_up - STATES.ca_rel)./CONSTANTS.tau_tr;
    ALGEBRAIC.i_up = CONSTANTS.i_up_max./(1.00000+CONSTANTS.k_up./STATES.cai);
    ALGEBRAIC.i_up_leak =  CONSTANTS.i_up_max.*(STATES.ca_up./CONSTANTS.ca_up_max);
    ALGEBRAIC.j_ca_csqn = RATES.dca_csqndt;
    ALGEBRAIC.j_ca_cmdn = RATES.dca_cmdndt;
    ALGEBRAIC.j_ca_trpn = RATES.dca_trpndt;
    ALGEBRAIC.qca_inf = 1.00000 - ((1.00000+((ALGEBRAIC.Fn./1.10000e-10) .^ 3.00000)) .^ -1.00000);
end

% Compute result of a piecewise function
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

function S1 = statesconv(S2)
if isstruct(S2)
    S1(:,1) = S2.V;
    S1(:,2) = S2.m;
    S1(:,3) = S2.h;
    S1(:,4) = S2.j;
    S1(:,5) = S2.d;
    S1(:,6)  = S2.f;
    S1(:,7)  = S2.xr;
    S1(:,8) = S2.xs;
    S1(:,9) = S2.nai;
    S1(:,10) = S2.cai;
    S1(:,11)  = S2.ki;
    S1(:,12) = S2.cli;
    S1(:,13) = S2.ca_up;
    S1(:,14)  = S2.ca_rel;
    S1(:,15) = S2.oa;
    S1(:,16) = S2.oi;
    S1(:,17)  = S2.qca;
    S1(:,18) = S2.ua;
    S1(:,19)  = S2.ui;
    S1(:,20) = S2.fca;
    S1(:,21) = S2.ca_cmdn;
    S1(:,22) = S2.ca_trpn;
    S1(:,23) = S2.ca_csqn;
    S1(:,24) = S2.u;
    S1(:,25)  = S2.v;
    S1(:,26)  = S2.w;
else
    S1.V = S2(:,1);
    S1.m = S2(:,2);
    S1.h = S2(:,3);
    S1.j = S2(:,4);
    S1.d = S2(:,5);
    S1.f = S2(:,6);
    S1.xr = S2(:,7);
    S1.xs = S2(:,8);
    S1.nai = S2(:,9);
    S1.cai = S2(:,10);
    S1.ki = S2(:,11);
    S1.cli = S2(:,12);
    S1.ca_up = S2(:,13);
    S1.ca_rel = S2(:,14);
    S1.oa = S2(:,15);
    S1.oi = S2(:,16);
    S1.qca= S2(:,17);
    S1.ua = S2(:,18);
    S1.ui = S2(:,19);
    S1.fca = S2(:,20);
    S1.ca_cmdn = S2(:,21);
    S1.ca_trpn = S2(:,22);
    S1.ca_csqn = S2(:,23);
    S1.u = S2(:,24);
    S1.v = S2(:,25);
    S1.w = S2(:,26);
end
end


function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (millisecond)');
    
    LEGEND_STATES.V = strpad('V in component membrane (millivolt)');
    LEGEND_STATES.nai = strpad('Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES.m = strpad('m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_STATES.h = strpad('h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_STATES.j = strpad('j in component fast_sodium_current_j_gate (dimensionless)');
    LEGEND_STATES.ki = strpad('K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES.oa = strpad('oa in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_STATES.oi = strpad('oi in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_STATES.ua = strpad('ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_STATES.ui = strpad('ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_STATES.xr = strpad('xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_STATES.xs = strpad('xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_STATES.cai = strpad('Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES.d = strpad('d in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_STATES.f = strpad('f in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_STATES.fca = strpad('f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_STATES.cli = strpad('Cl_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES.ca_rel = strpad('Ca_rel in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES.u = strpad('u in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_STATES.v = strpad('v in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_STATES.w = strpad('w in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_STATES.ca_cmdn = strpad('Ca_CMDN in component Ca_buffers (millimolar)');
    LEGEND_STATES.ca_trpn = strpad('Ca_TRPN in component Ca_buffers (millimolar)');
    LEGEND_STATES.ca_csqn = strpad('Ca_CSQN in component Ca_buffers (millimolar)');
    LEGEND_STATES.ca_up = strpad('Ca_up in component intracellular_ion_concentrations (millimolar)');
    
    LEGEND_RATES.dvdt = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES.dmdt = strpad('d/dt m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_RATES.dhdt = strpad('d/dt h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_RATES.djdt = strpad('d/dt j in component fast_sodium_current_j_gate (dimensionless)');
    LEGEND_RATES.doadt = strpad('d/dt oa in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_RATES.doidt = strpad('d/dt oi in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_RATES.duadt = strpad('d/dt ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_RATES.duidt = strpad('d/dt ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_RATES.dxrdt = strpad('d/dt xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_RATES.dxsdt = strpad('d/dt xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_RATES.dddt = strpad('d/dt d in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_RATES.dfdt = strpad('d/dt f in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_RATES.dfcadt = strpad('d/dt f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_RATES.dudt = strpad('d/dt u in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_RATES.djsrvdt = strpad('d/dt v in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_RATES.dwdt = strpad('d/dt w in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_RATES.dca_cmdndt = strpad('d/dt Ca_CMDN in component Ca_buffers (millimolar)');
    LEGEND_RATES.dca_trpndt = strpad('d/dt Ca_TRPN in component Ca_buffers (millimolar)');
    LEGEND_RATES.dca_csqndt = strpad('d/dt Ca_CSQN in component Ca_buffers (millimolar)');
    LEGEND_RATES.dnaidt = strpad('d/dt Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES.dkidt = strpad('d/dt K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES.dclidt = strpad('d/dt Cl_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES.dcaidt = strpad('d/dt Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES.dcaupdt = strpad('d/dt Ca_up in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES.dcareldt = strpad('d/dt Ca_rel in component intracellular_ion_concentrations (millimolar)');
    
    LEGEND_CONSTANTS.R = strpad('R in component membrane (joule_per_mole_kelvin)');
    LEGEND_CONSTANTS.T = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS.F = strpad('F in component membrane (coulomb_per_millimole)');
    LEGEND_CONSTANTS.Cm = strpad('Cm in component membrane (picoF)');
    LEGEND_CONSTANTS.istim = strpad('i_stim in component membrane (picoA_per_picoF)');
    LEGEND_CONSTANTS.gna = strpad('g_Na in component fast_sodium_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.nao = strpad('Na_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS.gk1 = strpad('g_K1 in component time_independent_potassium_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.ko = strpad('K_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS.gto = strpad('g_to in component transient_outward_K_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.gkr = strpad('g_Kr in component rapid_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.gks = strpad('g_Ks in component slow_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.gca = strpad('g_Ca in component sarcolemmal_Ca_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.tau_fca = strpad('tau_f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (millisecond)');
    LEGEND_CONSTANTS.gclca = strpad('g_Cl_Ca in component Ca_activated_Cl_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS.clo = strpad('Cl_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS.tau_qca = strpad('tau_qCa in component Ca_activated_Cl_current_q_Ca_gate (dimensionless)');
    LEGEND_CONSTANTS.kmnai = strpad('Km_Na_i in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS.kmko = strpad('Km_K_o in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS.inakmax = strpad('i_NaK_max in component sodium_potassium_pump (picoA_per_picoF)');
    LEGEND_CONSTANTS.sigma = strpad('sigma in component sodium_potassium_pump (dimensionless)');
    LEGEND_CONSTANTS.inacamax = strpad('I_NaCa_max in component Na_Ca_exchanger_current (picoA_per_picoF)');
    LEGEND_CONSTANTS.kmna = strpad('K_mNa in component Na_Ca_exchanger_current (millimolar)');
    LEGEND_CONSTANTS.kmca = strpad('K_mCa in component Na_Ca_exchanger_current (millimolar)');
    LEGEND_CONSTANTS.ksat = strpad('K_sat in component Na_Ca_exchanger_current (dimensionless)');
    LEGEND_CONSTANTS.cao = strpad('Ca_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS.g_b_na = strpad('g_B_Na in component background_currents (nanoS_per_picoF)');
    LEGEND_CONSTANTS.g_b_ca = strpad('g_B_Ca in component background_currents (nanoS_per_picoF)');
    LEGEND_CONSTANTS.i_p_ca_max = strpad('i_p_Ca_max in component Ca_pump_current (picoA_per_picoF)');
    LEGEND_CONSTANTS.k_rel = strpad('K_rel in component Ca_release_current_from_JSR (per_millisecond)');
    LEGEND_CONSTANTS.v_rel = strpad('V_rel in component Ca_release_current_from_JSR (micrometre_3)');
    LEGEND_CONSTANTS.tau_u = strpad('tau_u in component Ca_release_current_from_JSR_u_gate (millisecond)');
    LEGEND_CONSTANTS.tau_tr = strpad('tau_tr in component transfer_current_from_NSR_to_JSR (millisecond)');
    LEGEND_CONSTANTS.i_up_max = strpad('I_up_max in component Ca_uptake_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_CONSTANTS.ca_up_max = strpad('Ca_up_max in component Ca_leak_current_by_the_NSR (millimolar)');
    LEGEND_CONSTANTS.cmdn_max = strpad('CMDN_max in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS.trpn_max = strpad('TRPN_max in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS.csqn_max = strpad('CSQN_max in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS.k_up = strpad('K_up in component Ca_uptake_current_by_the_NSR (millimolar)');
    LEGEND_CONSTANTS.vi = strpad('V_i in component intracellular_ion_concentrations (micrometre_3)');
    LEGEND_CONSTANTS.v_rel = strpad('V_rel in component intracellular_ion_concentrations (micrometre_3)');
    LEGEND_CONSTANTS.v_up = strpad('V_up in component intracellular_ion_concentrations (micrometre_3)');
    
    LEGEND_ALGEBRAIC.ina = strpad('i_Na in component fast_sodium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ik1 = strpad('i_K1 in component time_independent_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ito = strpad('i_to in component transient_outward_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ikurd = strpad('i_Kur_d in component ultrarapid_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ikr = strpad('i_Kr in component rapid_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.iks = strpad('i_Ks in component slow_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ica = strpad('i_Ca in component sarcolemmal_Ca_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.iclca = strpad('i_Cl_Ca in component Ca_activated_Cl_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ipca = strpad('i_p_Ca in component Ca_pump_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.inak = strpad('i_NaK in component sodium_potassium_pump (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.inaca = strpad('i_NaCa in component Na_Ca_exchanger_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ibna = strpad('i_B_Na in component background_currents (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.ibca = strpad('i_B_Ca in component background_currents (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.Ena = strpad('E_Na in component fast_sodium_current (millivolt)');
    LEGEND_ALGEBRAIC.alpha_m = strpad('alpha_m in component fast_sodium_current_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_m = strpad('beta_m in component fast_sodium_current_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.alpha_h = strpad('alpha_h in component fast_sodium_current_h_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_h = strpad('beta_h in component fast_sodium_current_h_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.alpha_j = strpad('alpha_j in component fast_sodium_current_j_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_j = strpad('beta_j in component fast_sodium_current_j_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.Ek = strpad('E_K in component time_independent_potassium_current (millivolt)');
    LEGEND_ALGEBRAIC.alpha_oa = strpad('alpha_oa in component transient_outward_K_current_oa_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_oa = strpad('beta_oa in component transient_outward_K_current_oa_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_oa = strpad('tau_oa in component transient_outward_K_current_oa_gate (millisecond)');
    LEGEND_ALGEBRAIC.oa_inf = strpad('oa_infinity in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_ALGEBRAIC.alpha_oi = strpad('alpha_oi in component transient_outward_K_current_oi_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_oi = strpad('beta_oi in component transient_outward_K_current_oi_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_oi = strpad('tau_oi in component transient_outward_K_current_oi_gate (millisecond)');
    LEGEND_ALGEBRAIC.oi_inf = strpad('oi_infinity in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_ALGEBRAIC.gkurd = strpad('g_Kur_d in component ultrarapid_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_ALGEBRAIC.alpha_ua = strpad('alpha_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_ua = strpad('beta_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_ua = strpad('tau_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (millisecond)');
    LEGEND_ALGEBRAIC.ua_inf = strpad('ua_infinity in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_ALGEBRAIC.alpha_ui = strpad('alpha_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_ui = strpad('beta_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_ui = strpad('tau_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (millisecond)');
    LEGEND_ALGEBRAIC.ui_inf = strpad('ui_infinity in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_ALGEBRAIC.alpha_xr = strpad('alpha_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_xr = strpad('beta_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_xr = strpad('tau_xr in component rapid_delayed_rectifier_K_current_xr_gate (millisecond)');
    LEGEND_ALGEBRAIC.xr_inf = strpad('xr_infinity in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_ALGEBRAIC.alpha_xs = strpad('alpha_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.beta_xs = strpad('beta_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond)');
    LEGEND_ALGEBRAIC.tau_xs = strpad('tau_xs in component slow_delayed_rectifier_K_current_xs_gate (millisecond)');
    LEGEND_ALGEBRAIC.xs_inf = strpad('xs_infinity in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_ALGEBRAIC.dinf = strpad('d_infinity in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_ALGEBRAIC.tau_d = strpad('tau_d in component sarcolemmal_Ca_current_d_gate (millisecond)');
    LEGEND_ALGEBRAIC.f_inf = strpad('f_infinity in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_ALGEBRAIC.tau_f = strpad('tau_f in component sarcolemmal_Ca_current_f_gate (millisecond)');
    LEGEND_ALGEBRAIC.fca_inf = strpad('f_Ca_infinity in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_ALGEBRAIC.Ecl = strpad('E_Cl in component Ca_activated_Cl_current (millivolt)');
    LEGEND_ALGEBRAIC.Fn = strpad('Fn in component Ca_release_current_from_JSR (dimensionless)');
    LEGEND_ALGEBRAIC.qca_inf = strpad('q_Ca_infinity in component Ca_activated_Cl_current_q_Ca_gate (dimensionless)');
    LEGEND_ALGEBRAIC.fnak = strpad('f_NaK in component sodium_potassium_pump (dimensionless)');
    LEGEND_ALGEBRAIC.i_rel = strpad('i_rel in component Ca_release_current_from_JSR (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.Eca = strpad('E_Ca in component background_currents (millivolt)');
    LEGEND_ALGEBRAIC.u_inf = strpad('u_infinity in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_ALGEBRAIC.tau_v = strpad('tau_v in component Ca_release_current_from_JSR_v_gate (millisecond)');
    LEGEND_ALGEBRAIC.v_inf = strpad('v_infinity in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_ALGEBRAIC.tau_w = strpad('tau_w in component Ca_release_current_from_JSR_w_gate (millisecond)');
    LEGEND_ALGEBRAIC.w_inf = strpad('w_infinity in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_ALGEBRAIC.i_tr = strpad('i_tr in component transfer_current_from_NSR_to_JSR (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.i_up = strpad('i_up in component Ca_uptake_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.i_up_leak = strpad('i_up_leak in component Ca_leak_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_ALGEBRAIC.j_ca_cmdn = strpad('J_Ca_CMDN in component Ca_buffers (millimolar_per_millisecond)');
    LEGEND_ALGEBRAIC.j_ca_trpn = strpad('J_Ca_TRPN in component Ca_buffers (millimolar_per_millisecond)');
    LEGEND_ALGEBRAIC.j_ca_csqn = strpad('J_Ca_CSQN in component Ca_buffers (millimolar_per_millisecond)');

    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =68;
end
% There are a total of 25 entries in each of the rate and state variable arrays.
% There are a total of 42 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 1000];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-08, 'AbsTol', 1e-08, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, RATES', VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure(1);
    subplot(4,3,1); plot(VOI, STATES(:,1));  axis tight; xlim([0 300]); ylabel('AP');% V
    subplot(4,3,2); plot(VOI, ALGEBRAIC(:,15));  axis tight; xlim([0 300]); ylabel('Ina'); % Ina
    subplot(4,3,3); plot(VOI, ALGEBRAIC(:,47));   axis tight; xlim([0 300]); ylabel('Ica');% Ica
    subplot(4,3,4); plot(VOI, ALGEBRAIC(:,42));   axis tight; xlim([0 300]); ylabel('Ito');% Ito
    subplot(4,3,5); plot(VOI, ALGEBRAIC(:,45));  axis tight; xlim([0 300]); ylabel('Ikr');% Ikr
    subplot(4,3,6); plot(VOI, ALGEBRAIC(:,52));   axis tight; xlim([0 300]); ylabel('Inaca');% INaCa
    subplot(4,3,7); plot(VOI, ALGEBRAIC(:,44));   axis tight; xlim([0 300]); ylabel('Ikurd');% Ikurd
    subplot(4,3,8); plot(VOI, ALGEBRAIC(:,46));   axis tight; xlim([0 300]); ylabel('Iks');% Iks
    subplot(4,3,9); plot(VOI, ALGEBRAIC(:,51));   axis tight; xlim([0 300]); ylabel('Inak');% INak
    subplot(4,3,10); plot(VOI, ALGEBRAIC(:,49));  axis tight; xlim([0 300]);  ylabel('Iclca');% IClCa
    subplot(4,3,11); plot(VOI, ALGEBRAIC(:,35));   axis tight; xlim([0 300]); ylabel('Ik1');% Ik1
    subplot(4,3,12); plot(VOI, STATES(:,13)*1000);  ylabel('Cai');% Cai
    %     xlabel(LEGEND_VOI);
    %     l = legend(LEGEND_STATES);
    %     set(l,'Interpreter','none');
end


function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -83.53;
    CONSTANTS(:,1) = 8.3143;
    CONSTANTS(:,2) = 310.0;
    CONSTANTS(:,3) = 96.4867;
    CONSTANTS(:,4) = 1; % Cm = 100pf
    CONSTANTS(:,5) = -80.0; %stim_amplitude in component membrane (picoA_per_picoF)
    CONSTANTS(:,6) = 8.8;  %gna = 7.8
    STATES(:,2) = 11.75;
    CONSTANTS(:,7) = 140.0;
    STATES(:,3) = 0.001972;
    STATES(:,4) = 0.9791;
    STATES(:,5) = 0.9869;
    CONSTANTS(:,8) = 0.15;
    CONSTANTS(:,9) = 5.4; % ko = 5.4
    STATES(:,6) = 138.4;
    CONSTANTS(:,10) = 0.19824;
    STATES(:,7) = 0.07164;
    STATES(:,8) = 0.9980;
    STATES(:,9) = 0.05869;
    STATES(:,10) = 0.9987;
    CONSTANTS(:,11) = 0.0105;  %gkr = 0.06984
    STATES(:,11) = 0.0000007433;
    CONSTANTS(:,12) = 0.0561;
    STATES(:,12) = 0.01791;
    CONSTANTS(:,13) = 0.24; %gca = 0.24
    STATES(:,13) = 0.0001024; %cai= 0.0001024
    STATES(:,14) = 0.000004757;
    STATES(:,15) = 0.9999;
    STATES(:,16) = 0.7484;
    CONSTANTS(:,14) = 2.0; %tau_fca=2ms
    CONSTANTS(:,15) = 0.3; %gclca = 0.3;
    STATES(:,17) = 29.26;
    CONSTANTS(:,16) = 132.0;
    CONSTANTS(:,17) = 0.0;
    CONSTANTS(:,18) = 10.0;
    CONSTANTS(:,19) = 1.5;
    CONSTANTS(:,20) = 0.6;
    CONSTANTS(:,21) = 1600.0;
    CONSTANTS(:,22) = 87.5;
    CONSTANTS(:,23) = 1.38;
    CONSTANTS(:,24) = 0.1;
    CONSTANTS(:,25) = 1.8; %cao=1.8
    CONSTANTS(:,26) = 0.000674;
    CONSTANTS(:,27) = 0.00113;
    CONSTANTS(:,28) = 0.275;
    CONSTANTS(:,29) = 30.0;
    CONSTANTS(:,30) = 96.48;
    STATES(:,18) = 1.502;
    STATES(:,19) = 0.0;
    STATES(:,20) = 1.0;
    STATES(:,21) = 0.9993;
    CONSTANTS(:,31) = 180.0;
    STATES(:,22) = 1.502;
    CONSTANTS(:,32) = 0.005;
    CONSTANTS(:,33) = 0.00092;
    CONSTANTS(:,34) = 15.0;   %ca_up_max = 15.0;
    CONSTANTS(:,35) = 0.045;
    CONSTANTS(:,36) = 0.35;
    CONSTANTS(:,37) = 10.0;
    STATES(:,23) = 0.001856;
    STATES(:,24) = 0.007022;
    STATES(:,25) = 6.432;
    CONSTANTS(:,38) = 13668.0;
    CONSTANTS(:,39) = 96.48;
    CONSTANTS(:,40) = 1109.52;
    CONSTANTS(:,41) =  (1.00000./7.00000).*((exp((CONSTANTS(:,7)./67.3000))) - 1.00000);
    CONSTANTS(:,42) = 8.00000;
    CONSTANTS(:,43) = 50; %stim_start in component membrane (millisecond)
    CONSTANTS(:,44) = 1000; %stim_period in component membrane (millisecond)
    CONSTANTS(:,45) = 1; %stim_duration in component membrane (millisecond)
    if (isempty(STATES))
        warning('Initial values for states not set'); 
    end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
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
    
    ALGEBRAIC(:,13) = 0.290000+ 0.800000.*((1.00000+(exp(((STATES(:,13) - 0.000120000)./6.00000e-05)))) .^ -1.00000);
    RATES(:,16) = (ALGEBRAIC(:,13) - STATES(:,16))./CONSTANTS(:,14);
    ALGEBRAIC(:,2) =  0.320000.*((STATES(:,1)+47.1300)./(1.00000 - (exp(( -0.100000.*(STATES(:,1)+47.1300))))));
    ALGEBRAIC(:,16) =  0.0800000.*(exp((STATES(:,1)./-11.0000)));
    RATES(:,3) =  ALGEBRAIC(:,2).*(1.00000 - STATES(:,3)) -  ALGEBRAIC(:,16).*STATES(:,3);
    ALGEBRAIC(:,3) = piecewise({STATES(:,1)<-40.0000,  0.135000.*(exp(((STATES(:,1)+80.0000)./-6.80000))) }, 0.00000);
    ALGEBRAIC(:,17) = piecewise({STATES(:,1)<-40.0000,  3.56000.*(exp(( 0.0790000.*STATES(:,1))))+ 310000..*(exp(( 0.350000.*STATES(:,1)))) }, 1.00000./( 0.130000.*(1.00000+(exp(((STATES(:,1)+10.6600)./-11.1000))))));
    RATES(:,4) =  ALGEBRAIC(:,3).*(1.00000 - STATES(:,4)) -  ALGEBRAIC(:,17).*STATES(:,4);
    ALGEBRAIC(:,4) = piecewise({STATES(:,1)<-40.0000,  (( -127140..*(exp(( 0.244400.*STATES(:,1)))) -  3.47400e-05.*(exp(( -0.0439100.*STATES(:,1)))))./(1.00000+(exp(( 0.311000.*(STATES(:,1)+79.2300)))))).*(STATES(:,1)+37.7800) }, 0.00000);
    ALGEBRAIC(:,18) = piecewise({STATES(:,1)<-40.0000, ( 0.121200.*(exp(( -0.0105200.*STATES(:,1)))))./(1.00000+(exp(( -0.137800.*(STATES(:,1)+40.1400))))) }, ( 0.300000.*(exp(( -2.53500e-07.*STATES(:,1)))))./(1.00000+(exp(( -0.100000.*(STATES(:,1)+32.0000))))));
    RATES(:,5) =  ALGEBRAIC(:,4).*(1.00000 - STATES(:,5)) -  ALGEBRAIC(:,18).*STATES(:,5);
    
    
    ALGEBRAIC(:,11) = (1.00000+(exp(((STATES(:,1)+10.0000)./-6.00000)))) .^ -1.00000;
    ALGEBRAIC(:,25) = (1.00000 - (exp(((STATES(:,1)+10.0000)./-6.24000))))./( 0.0350000.*(STATES(:,1)...
        +10.0000).*(1.00000+(exp(((STATES(:,1)+10.0000)./-6.24000)))));
    RATES(:,14) = (ALGEBRAIC(:,11) - STATES(:,14))./ALGEBRAIC(:,25);
    
    
    ALGEBRAIC(:,12) = (1.00000+(exp(((STATES(:,1)+24.6000)./6.20000)))) .^ -1.00000;
    ALGEBRAIC(:,26) =  400.000.*((1.00000+ 4.50000.*(exp(( -0.000700000.*((STATES(:,1) - 9.00000) .^ 2.00000))))) .^ -1.00000);
    RATES(:,15) = (ALGEBRAIC(:,12) - STATES(:,15))./ALGEBRAIC(:,26);
    ALGEBRAIC(:,14) = (6.00000 -  6.00000.*(exp(((STATES(:,1) - 7.90000)./-5.00000))))./( (1.00000+ 0.300000.*(exp(((STATES(:,1) - 7.90000)./-5.00000)))).*(STATES(:,1) - 7.90000));
    ALGEBRAIC(:,27) = 1.00000 - ((1.00000+(exp(((STATES(:,1) - 40.0000)./-17.0000)))) .^ -1.00000);
    RATES(:,21) = (ALGEBRAIC(:,27) - STATES(:,21))./ALGEBRAIC(:,14);
    ALGEBRAIC(:,5) =  0.650000.*(((exp(((STATES(:,1)+18.0000)./-8.50000)))+(exp(((STATES(:,1) - 16.0000)./-59.0000)))) .^ -1.00000);
    ALGEBRAIC(:,19) =  1.20000.*((2.20000+(exp(((STATES(:,1)+75.0000)./18.0000)))) .^ -1.00000);
    ALGEBRAIC(:,29) = (ALGEBRAIC(:,5)+ALGEBRAIC(:,19)) .^ -1.00000;
    ALGEBRAIC(:,36) = (1.00000+(exp(((STATES(:,1)+0.500000)./-10.5000)))) .^ ( - (1.00000./3.00000));
    RATES(:,7) = (ALGEBRAIC(:,36) - STATES(:,7))./ALGEBRAIC(:,29);
    ALGEBRAIC(:,6) = (6.20000+(exp(((STATES(:,1)+105.200)./9.85000)))) .^ -1.00000;
    ALGEBRAIC(:,20) = (7.54000+(exp(((STATES(:,1) - 8.89000)./-12.8700)))) .^ -1.00000;
    ALGEBRAIC(:,30) = (ALGEBRAIC(:,6)+ALGEBRAIC(:,20)) .^ -1.00000;
    ALGEBRAIC(:,37) = (1.00000+(exp(((STATES(:,1)+43.3770)./6.45000)))) .^ -1.00000;
    RATES(:,8) = (ALGEBRAIC(:,37) - STATES(:,8))./ALGEBRAIC(:,30);
    ALGEBRAIC(:,7) =  1.47000.*(((exp(((STATES(:,1)+33.2000)./-30.6300)))+(exp(((STATES(:,1) - 27.6000)./-30.6500)))) .^ -1.00000);
    ALGEBRAIC(:,21) =  0.420000.*(((exp(((STATES(:,1)+26.6400)./2.49000)))+(exp(((STATES(:,1)+44.4100)./20.3600)))) .^ -1.00000);
    ALGEBRAIC(:,31) = (ALGEBRAIC(:,7)+ALGEBRAIC(:,21)) .^ -1.00000;
    ALGEBRAIC(:,38) = (1.00000+(exp(((STATES(:,1)+2.81000)./-9.49000)))) .^ ( - (1.00000./3.00000));
    RATES(:,9) = (ALGEBRAIC(:,38) - STATES(:,9))./ALGEBRAIC(:,31);
    ALGEBRAIC(:,8) = (21.0000+(exp(((STATES(:,1) - 185.000)./-28.0000)))) .^ -1.00000;
    ALGEBRAIC(:,22) = exp(((STATES(:,1) - 158.000)./16.0000));
    ALGEBRAIC(:,32) = (ALGEBRAIC(:,8)+ALGEBRAIC(:,22)) .^ -1.00000;
    ALGEBRAIC(:,39) = (1.00000+(exp(((STATES(:,1) - 99.4500)./27.4800)))) .^ -1.00000;
    RATES(:,10) = (ALGEBRAIC(:,39) - STATES(:,10))./ALGEBRAIC(:,32);
    ALGEBRAIC(:,9) =  0.0400000.*((STATES(:,1) - 248.000)./(1.00000 - (exp(((STATES(:,1) - 248.000)./-28.0000)))));
    ALGEBRAIC(:,23) =  0.0280000.*((STATES(:,1)+163.000)./((exp(((STATES(:,1)+163.000)./21.0000))) - 1.00000));
    ALGEBRAIC(:,33) = (ALGEBRAIC(:,9)+ALGEBRAIC(:,23)) .^ -1.00000;
    ALGEBRAIC(:,40) = (1.00000+(exp(((STATES(:,1)+7.65400)./-5.37700)))) .^ -1.00000;
    RATES(:,11) = (ALGEBRAIC(:,40) - STATES(:,11))./ALGEBRAIC(:,33);
    ALGEBRAIC(:,10) =  1.00000e-05.*((STATES(:,1)+28.5000)./(1.00000 - (exp(((STATES(:,1)+28.5000)./-115.000)))));
    ALGEBRAIC(:,24) =  0.000230000.*((STATES(:,1)+28.5000)./((exp(((STATES(:,1)+28.5000)./3.30000))) - 1.00000));
    ALGEBRAIC(:,34) = (ALGEBRAIC(:,10)+ALGEBRAIC(:,24)) .^ -1.00000;
    ALGEBRAIC(:,41) = (1.00000+(exp(((STATES(:,1) - 13.0000)./-12.0000)))) .^ -0.500000;
    RATES(:,12) = (ALGEBRAIC(:,41) - STATES(:,12))./ALGEBRAIC(:,34);
    
    
    ALGEBRAIC(:,48) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( -1.00000.*CONSTANTS(:,3))).*(log((CONSTANTS(:,16)./STATES(:,17))));
    ALGEBRAIC(:,49) =  CONSTANTS(:,15).*CONSTANTS(:,17).*(STATES(:,1) - ALGEBRAIC(:,48)); %iclca
    
    
    RATES(:,17) = ALGEBRAIC(:,49)./( CONSTANTS(:,38).*CONSTANTS(:,3));
    ALGEBRAIC(:,28) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*(log((CONSTANTS(:,9)./STATES(:,6))));
    ALGEBRAIC(:,35) = ( CONSTANTS(:,8).*(STATES(:,1) - ALGEBRAIC(:,28)))./(1.00000+(exp(( 0.0700000.*(STATES(:,1)+80.0000)))));
    ALGEBRAIC(:,42) =  CONSTANTS(:,10).*(STATES(:,7) .^ 3.00000).*STATES(:,8).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,43) = 0.00855000+0.0779000./(1.00000+(exp(((STATES(:,1)+11.0000)./-16.0000))));
    ALGEBRAIC(:,44) =  ALGEBRAIC(:,43).*(STATES(:,9) .^ 3.00000).*STATES(:,10).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,45) =  CONSTANTS(:,11).*STATES(:,11).*(0.0700000+0.580000./(1.00000+(exp(((STATES(:,1)+...
        15.0000)./22.4000))))).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,46) =  CONSTANTS(:,12).*(STATES(:,12) .^ 2.00000).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,50) = (1.00000+ 0.124500.*(exp(( -0.100000.*(( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))))+ 0.0365000.*CONSTANTS(:,41).*(exp(( - (( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))))) .^ -1.00000;
    ALGEBRAIC(:,51) =  CONSTANTS(:,20).*ALGEBRAIC(:,50).*(1.00000./(1.00000+((CONSTANTS(:,18)./STATES(:,2)) .^ 1.50000))).*(CONSTANTS(:,9)./(CONSTANTS(:,9)+CONSTANTS(:,19)));
    RATES(:,6) = ( 2.00000.*ALGEBRAIC(:,51) - (ALGEBRAIC(:,35)+ALGEBRAIC(:,42)+ALGEBRAIC(:,44)+ALGEBRAIC(:,45)+ALGEBRAIC(:,46)))./( CONSTANTS(:,38).*CONSTANTS(:,3));
    ALGEBRAIC(:,1) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*(log((CONSTANTS(:,7)./STATES(:,2))));
    ALGEBRAIC(:,15) =  CONSTANTS(:,6).*(STATES(:,3) .^ 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,52) = ( CONSTANTS(:,21).*( (exp((( 0.350000.*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*(STATES(:,2) .^ 3.00000).*CONSTANTS(:,25) -  (exp((( -0.650000.*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*(CONSTANTS(:,7) .^ 3.00000).*STATES(:,13)))./( ((CONSTANTS(:,22) .^ 3.00000)+(CONSTANTS(:,7) .^ 3.00000)).*(CONSTANTS(:,23)+CONSTANTS(:,25)).*(1.00000+ CONSTANTS(:,24).*(exp((( -0.650000.*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))))));
    ALGEBRAIC(:,53) =  CONSTANTS(:,26).*(STATES(:,1) - ALGEBRAIC(:,1));
    RATES(:,2) = ( -3.00000.*ALGEBRAIC(:,51) - ( 3.00000.*ALGEBRAIC(:,52)+ALGEBRAIC(:,53)+ALGEBRAIC(:,15)))./( CONSTANTS(:,38).*CONSTANTS(:,3));
    ALGEBRAIC(:,47) =  CONSTANTS(:,13).*STATES(:,14).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - 65.0000);
    ALGEBRAIC(:,56) =  CONSTANTS(:,28).*(STATES(:,13)./(0.000500000+STATES(:,13)));
    ALGEBRAIC(:,54) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*(log((CONSTANTS(:,25)./STATES(:,13))));
    ALGEBRAIC(:,55) =  CONSTANTS(:,27).*(STATES(:,1) - ALGEBRAIC(:,54));
    
    STIM = piecewise({VOI -  (floor((VOI./CONSTANTS(:,44)))).*CONSTANTS(:,44)>=CONSTANTS(:,43)& VOI -  ...
        (floor((VOI./CONSTANTS(:,44)))).*CONSTANTS(:,44)<=CONSTANTS(:,43)+CONSTANTS(:,45),  CONSTANTS(:,5)}, 0.00000);
    
    RATES(:,1) =  - (ALGEBRAIC(:,15)+ALGEBRAIC(:,35)+ALGEBRAIC(:,42)+ALGEBRAIC(:,44)+ALGEBRAIC(:,45)+ALGEBRAIC(:,46)+ALGEBRAIC(:,47)...
        +ALGEBRAIC(:,49)+ALGEBRAIC(:,56)+ALGEBRAIC(:,52)+ALGEBRAIC(:,51)+ALGEBRAIC(:,53)+ALGEBRAIC(:,55)+STIM)./CONSTANTS(:,4);
    
    RATES(:,25) =  0.480000.*STATES(:,18).*(1.00000 - STATES(:,25)./CONSTANTS(:,37)) -  0.400000.*(STATES(:,25)./CONSTANTS(:,37));
    ALGEBRAIC(:,57) =  CONSTANTS(:,29).*(STATES(:,19) .^ 2.00000).*STATES(:,20).*STATES(:,21).*(STATES(:,18) - STATES(:,13));
    ALGEBRAIC(:,58) =  1.00000e-12.*CONSTANTS(:,30).*ALGEBRAIC(:,57) -  5.00000e-13.*( (1.00000./( 2.00000.*CONSTANTS(:,3))).*ALGEBRAIC(:,47) -  (1.00000./( 5.00000.*CONSTANTS(:,3))).*ALGEBRAIC(:,52));
    ALGEBRAIC(:,61) = (1.00000+(exp(((ALGEBRAIC(:,58) - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000;
    RATES(:,19) = (ALGEBRAIC(:,61) - STATES(:,19))./CONSTANTS(:,42);
    ALGEBRAIC(:,62) = 1.91000+ 2.09000.*((1.00000+(exp(((ALGEBRAIC(:,58) - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000);
    ALGEBRAIC(:,64) = 1.00000 - ((1.00000+(exp(((ALGEBRAIC(:,58) - 6.83500e-14)./-1.36700e-15)))) .^ -1.00000);
    RATES(:,20) = (ALGEBRAIC(:,64) - STATES(:,20))./ALGEBRAIC(:,62);
    RATES(:,23) =  200.000.*STATES(:,13).*(1.00000 - STATES(:,23)./CONSTANTS(:,35)) -  0.476000.*(STATES(:,23)./CONSTANTS(:,35));
    ALGEBRAIC(:,59) = (STATES(:,22) - STATES(:,18))./CONSTANTS(:,31);
    ALGEBRAIC(:,63) = CONSTANTS(:,32)./(1.00000+CONSTANTS(:,33)./STATES(:,13));
    ALGEBRAIC(:,66) =  CONSTANTS(:,32).*(STATES(:,22)./CONSTANTS(:,34));
    RATES(:,22) = ALGEBRAIC(:,63) - (ALGEBRAIC(:,66)+ ALGEBRAIC(:,59).*(CONSTANTS(:,39)./CONSTANTS(:,40)));
    ALGEBRAIC(:,65) = RATES(:,25);
    RATES(:,18) = ALGEBRAIC(:,59) - (ALGEBRAIC(:,57)+ 31.0000.*ALGEBRAIC(:,65));
    RATES(:,24) =  78.4000.*STATES(:,13).*(1.00000 - STATES(:,24)./CONSTANTS(:,36)) -  0.392000.*(STATES(:,24)./CONSTANTS(:,36));
    ALGEBRAIC(:,67) = RATES(:,23);
    ALGEBRAIC(:,68) = RATES(:,24);
    RATES(:,13) = (( 2.00000.*ALGEBRAIC(:,52) - (ALGEBRAIC(:,56)+ALGEBRAIC(:,47)+ALGEBRAIC(:,55)))./( 2.00000.*CONSTANTS(:,38).*CONSTANTS(:,3))+( CONSTANTS(:,40).*(ALGEBRAIC(:,66) - ALGEBRAIC(:,63))+ ALGEBRAIC(:,57).*CONSTANTS(:,39))./CONSTANTS(:,38)) - ( CONSTANTS(:,36).*ALGEBRAIC(:,68)+ CONSTANTS(:,35).*ALGEBRAIC(:,67));
    RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, RATES, VOI)
    ALGEBRAIC(:,13) = 0.290000+ 0.800000.*((1.00000+(exp(((STATES(:,13) - 0.000120000)./6.00000e-05)))) .^ -1.00000);
    ALGEBRAIC(:,2) =  0.320000.*((STATES(:,1)+47.1300)./(1.00000 - (exp(( -0.100000.*(STATES(:,1)+47.1300))))));
    ALGEBRAIC(:,16) =  0.0800000.*(exp((STATES(:,1)./-11.0000)));
    ALGEBRAIC(:,3) = piecewise({STATES(:,1)<-40.0000,  0.135000.*(exp(((STATES(:,1)+80.0000)./-6.80000))) }, 0.00000);
    ALGEBRAIC(:,17) = piecewise({STATES(:,1)<-40.0000,  3.56000.*(exp(( 0.0790000.*STATES(:,1))))+ 310000..*(exp(( 0.350000.*STATES(:,1)))) }, 1.00000./( 0.130000.*(1.00000+(exp(((STATES(:,1)+10.6600)./-11.1000))))));
    ALGEBRAIC(:,4) = piecewise({STATES(:,1)<-40.0000,  (( -127140..*(exp(( 0.244400.*STATES(:,1)))) -  3.47400e-05.*(exp(( -0.0439100.*STATES(:,1)))))./(1.00000+(exp(( 0.311000.*(STATES(:,1)+79.2300)))))).*(STATES(:,1)+37.7800) }, 0.00000);
    ALGEBRAIC(:,18) = piecewise({STATES(:,1)<-40.0000, ( 0.121200.*(exp(( -0.0105200.*STATES(:,1)))))./(1.00000+(exp(( -0.137800.*(STATES(:,1)+40.1400))))) }, ( 0.300000.*(exp(( -2.53500e-07.*STATES(:,1)))))./(1.00000+(exp(( -0.100000.*(STATES(:,1)+32.0000))))));
    ALGEBRAIC(:,11) = (1.00000+(exp(((STATES(:,1)+10.0000)./-6.00000)))) .^ -1.00000;
    ALGEBRAIC(:,25) = (1.00000 - (exp(((STATES(:,1)+10.0000)./-6.24000))))./( 0.0350000.*(STATES(:,1)+10.0000).*(1.00000+(exp(((STATES(:,1)+10.0000)./-6.24000)))));
    ALGEBRAIC(:,12) = (1.00000+(exp(((STATES(:,1)+24.6000)./6.20000)))) .^ -1.00000;
    ALGEBRAIC(:,26) =  400.000.*((1.00000+ 4.50000.*(exp(( -0.000700000.*((STATES(:,1) - 9.00000) .^ 2.00000))))) .^ -1.00000);
    ALGEBRAIC(:,14) = (6.00000 -  6.00000.*(exp(((STATES(:,1) - 7.90000)./-5.00000))))./( (1.00000+ 0.300000.*(exp(((STATES(:,1) - 7.90000)./-5.00000)))).*(STATES(:,1) - 7.90000));
    ALGEBRAIC(:,27) = 1.00000 - ((1.00000+(exp(((STATES(:,1) - 40.0000)./-17.0000)))) .^ -1.00000);
    ALGEBRAIC(:,5) =  0.650000.*(((exp(((STATES(:,1)+18.0000)./-8.50000)))+(exp(((STATES(:,1) - 16.0000)./-59.0000)))) .^ -1.00000);
    ALGEBRAIC(:,19) =  1.20000.*((2.20000+(exp(((STATES(:,1)+75.0000)./18.0000)))) .^ -1.00000);
    ALGEBRAIC(:,29) = (ALGEBRAIC(:,5)+ALGEBRAIC(:,19)) .^ -1.00000;
    ALGEBRAIC(:,36) = (1.00000+(exp(((STATES(:,1)+0.500000)./-10.5000)))) .^ ( - (1.00000./3.00000));
    ALGEBRAIC(:,6) = (6.20000+(exp(((STATES(:,1)+105.200)./9.85000)))) .^ -1.00000;
    ALGEBRAIC(:,20) = (7.54000+(exp(((STATES(:,1) - 8.89000)./-12.8700)))) .^ -1.00000;
    ALGEBRAIC(:,30) = (ALGEBRAIC(:,6)+ALGEBRAIC(:,20)) .^ -1.00000;
    ALGEBRAIC(:,37) = (1.00000+(exp(((STATES(:,1)+43.3770)./6.45000)))) .^ -1.00000;
    ALGEBRAIC(:,7) =  1.47000.*(((exp(((STATES(:,1)+33.2000)./-30.6300)))+(exp(((STATES(:,1) - 27.6000)./-30.6500)))) .^ -1.00000);
    ALGEBRAIC(:,21) =  0.420000.*(((exp(((STATES(:,1)+26.6400)./2.49000)))+(exp(((STATES(:,1)+44.4100)./20.3600)))) .^ -1.00000);
    ALGEBRAIC(:,31) = (ALGEBRAIC(:,7)+ALGEBRAIC(:,21)) .^ -1.00000;
    ALGEBRAIC(:,38) = (1.00000+(exp(((STATES(:,1)+2.81000)./-9.49000)))) .^ ( - (1.00000./3.00000));
    ALGEBRAIC(:,8) = (21.0000+(exp(((STATES(:,1) - 185.000)./-28.0000)))) .^ -1.00000;
    ALGEBRAIC(:,22) = exp(((STATES(:,1) - 158.000)./16.0000));
    ALGEBRAIC(:,32) = (ALGEBRAIC(:,8)+ALGEBRAIC(:,22)) .^ -1.00000;
    ALGEBRAIC(:,39) = (1.00000+(exp(((STATES(:,1) - 99.4500)./27.4800)))) .^ -1.00000;
    ALGEBRAIC(:,9) =  0.0400000.*((STATES(:,1) - 248.000)./(1.00000 - (exp(((STATES(:,1) - 248.000)./-28.0000)))));
    ALGEBRAIC(:,23) =  0.0280000.*((STATES(:,1)+163.000)./((exp(((STATES(:,1)+163.000)./21.0000))) - 1.00000));
    ALGEBRAIC(:,33) = (ALGEBRAIC(:,9)+ALGEBRAIC(:,23)) .^ -1.00000;
    ALGEBRAIC(:,40) = (1.00000+(exp(((STATES(:,1)+7.65400)./-5.37700)))) .^ -1.00000;
    ALGEBRAIC(:,10) =  1.00000e-05.*((STATES(:,1)+28.5000)./(1.00000 - (exp(((STATES(:,1)+28.5000)./-115.000)))));
    ALGEBRAIC(:,24) =  0.000230000.*((STATES(:,1)+28.5000)./((exp(((STATES(:,1)+28.5000)./3.30000))) - 1.00000));
    ALGEBRAIC(:,34) = (ALGEBRAIC(:,10)+ALGEBRAIC(:,24)) .^ -1.00000;
    ALGEBRAIC(:,41) = (1.00000+(exp(((STATES(:,1) - 13.0000)./-12.0000)))) .^ -0.500000;
    ALGEBRAIC(:,48) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( -1.00000.*CONSTANTS(:,3))).*(log((CONSTANTS(:,16)./STATES(:,17))));
    ALGEBRAIC(:,49) =  CONSTANTS(:,15).*CONSTANTS(:,17).*(STATES(:,1) - ALGEBRAIC(:,48));
    ALGEBRAIC(:,28) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*(log((CONSTANTS(:,9)./STATES(:,6))));
    ALGEBRAIC(:,35) = ( CONSTANTS(:,8).*(STATES(:,1) - ALGEBRAIC(:,28)))./(1.00000+(exp(( 0.0700000.*(STATES(:,1)+80.0000)))));
    ALGEBRAIC(:,42) =  CONSTANTS(:,10).*(STATES(:,7) .^ 3.00000).*STATES(:,8).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,43) = 0.00855000+0.0779000./(1.00000+(exp(((STATES(:,1)+11.0000)./-16.0000))));
    ALGEBRAIC(:,44) =  ALGEBRAIC(:,43).*(STATES(:,9) .^ 3.00000).*STATES(:,10).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,45) =  CONSTANTS(:,11).*STATES(:,11).*(0.0700000+0.580000./(1.00000+(exp(((STATES(:,1)+15.0000)./22.4000))))).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,46) =  CONSTANTS(:,12).*(STATES(:,12) .^ 2.00000).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,50) = (1.00000+ 0.124500.*(exp(( -0.100000.*(( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))))+ 0.0365000.*CONSTANTS(:,41).*(exp(( - (( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))))) .^ -1.00000;
    ALGEBRAIC(:,51) =  CONSTANTS(:,20).*ALGEBRAIC(:,50).*(1.00000./(1.00000+((CONSTANTS(:,18)./STATES(:,2)) .^ 1.50000))).*(CONSTANTS(:,9)./(CONSTANTS(:,9)+CONSTANTS(:,19)));
    ALGEBRAIC(:,1) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*(log((CONSTANTS(:,7)./STATES(:,2))));
    ALGEBRAIC(:,15) =  CONSTANTS(:,6).*(STATES(:,3) .^ 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,52) = ( CONSTANTS(:,21).*( (exp((( 0.350000.*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*(STATES(:,2) .^ 3.00000).*CONSTANTS(:,25) -  (exp((( -0.650000.*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*(CONSTANTS(:,7) .^ 3.00000).*STATES(:,13)))./( ((CONSTANTS(:,22) .^ 3.00000)+(CONSTANTS(:,7) .^ 3.00000)).*(CONSTANTS(:,23)+CONSTANTS(:,25)).*(1.00000+ CONSTANTS(:,24).*(exp((( -0.650000.*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))))));
    ALGEBRAIC(:,53) =  CONSTANTS(:,26).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,47) =  CONSTANTS(:,13).*STATES(:,14).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - 65.0000);
    ALGEBRAIC(:,56) =  CONSTANTS(:,28).*(STATES(:,13)./(0.000500000+STATES(:,13)));
    ALGEBRAIC(:,54) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*(log((CONSTANTS(:,25)./STATES(:,13))));
    ALGEBRAIC(:,55) =  CONSTANTS(:,27).*(STATES(:,1) - ALGEBRAIC(:,54));
    ALGEBRAIC(:,57) =  CONSTANTS(:,29).*(STATES(:,19) .^ 2.00000).*STATES(:,20).*STATES(:,21).*(STATES(:,18) - STATES(:,13));
    ALGEBRAIC(:,58) =  1.00000e-12.*CONSTANTS(:,30).*ALGEBRAIC(:,57) -  5.00000e-13.*( (1.00000./( 2.00000.*CONSTANTS(:,3))).*ALGEBRAIC(:,47) -  (1.00000./( 5.00000.*CONSTANTS(:,3))).*ALGEBRAIC(:,52));
    ALGEBRAIC(:,61) = (1.00000+(exp(((ALGEBRAIC(:,58) - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000;
    ALGEBRAIC(:,62) = 1.91000+ 2.09000.*((1.00000+(exp(((ALGEBRAIC(:,58) - 3.41750e-13)./-1.36700e-15)))) .^ -1.00000);
    ALGEBRAIC(:,64) = 1.00000 - ((1.00000+(exp(((ALGEBRAIC(:,58) - 6.83500e-14)./-1.36700e-15)))) .^ -1.00000);
    ALGEBRAIC(:,59) = (STATES(:,22) - STATES(:,18))./CONSTANTS(:,31);
    ALGEBRAIC(:,63) = CONSTANTS(:,32)./(1.00000+CONSTANTS(:,33)./STATES(:,13));
    ALGEBRAIC(:,66) =  CONSTANTS(:,32).*(STATES(:,22)./CONSTANTS(:,34));
    ALGEBRAIC(:,65) = RATES(:,25);
    ALGEBRAIC(:,67) = RATES(:,23);
    ALGEBRAIC(:,68) = RATES(:,24);
    ALGEBRAIC(:,60) = 1.00000 - ((1.00000+((ALGEBRAIC(:,58)./1.10000e-10) .^ 3.00000)) .^ -1.00000);
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


function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (millisecond)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component membrane (joule_per_mole_kelvin)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component membrane (coulomb_per_millimole)');
    LEGEND_CONSTANTS(:,4) = strpad('Cm in component membrane (picoF)');
    LEGEND_CONSTANTS(:,5) = strpad('i_stim in component membrane (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,15) = strpad('i_Na in component fast_sodium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,35) = strpad('i_K1 in component time_independent_potassium_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,42) = strpad('i_to in component transient_outward_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,44) = strpad('i_Kur_d in component ultrarapid_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,45) = strpad('i_Kr in component rapid_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,46) = strpad('i_Ks in component slow_delayed_rectifier_K_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,47) = strpad('i_Ca in component sarcolemmal_Ca_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,49) = strpad('i_Cl_Ca in component Ca_activated_Cl_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,56) = strpad('i_p_Ca in component Ca_pump_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,51) = strpad('i_NaK in component sodium_potassium_pump (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,52) = strpad('i_NaCa in component Na_Ca_exchanger_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,53) = strpad('i_B_Na in component background_currents (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,55) = strpad('i_B_Ca in component background_currents (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,1) = strpad('E_Na in component fast_sodium_current (millivolt)');
    LEGEND_CONSTANTS(:,6) = strpad('g_Na in component fast_sodium_current (nanoS_per_picoF)');
    LEGEND_STATES(:,2) = strpad('Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,7) = strpad('Na_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_STATES(:,3) = strpad('m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_STATES(:,4) = strpad('h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_STATES(:,5) = strpad('j in component fast_sodium_current_j_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('alpha_m in component fast_sodium_current_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,16) = strpad('beta_m in component fast_sodium_current_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,3) = strpad('alpha_h in component fast_sodium_current_h_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,17) = strpad('beta_h in component fast_sodium_current_h_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,4) = strpad('alpha_j in component fast_sodium_current_j_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,18) = strpad('beta_j in component fast_sodium_current_j_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,28) = strpad('E_K in component time_independent_potassium_current (millivolt)');
    LEGEND_CONSTANTS(:,8) = strpad('g_K1 in component time_independent_potassium_current (nanoS_per_picoF)');
    LEGEND_CONSTANTS(:,9) = strpad('K_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_STATES(:,6) = strpad('K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,10) = strpad('g_to in component transient_outward_K_current (nanoS_per_picoF)');
    LEGEND_STATES(:,7) = strpad('oa in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_STATES(:,8) = strpad('oi in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('alpha_oa in component transient_outward_K_current_oa_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,19) = strpad('beta_oa in component transient_outward_K_current_oa_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,29) = strpad('tau_oa in component transient_outward_K_current_oa_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,36) = strpad('oa_infinity in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('alpha_oi in component transient_outward_K_current_oi_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,20) = strpad('beta_oi in component transient_outward_K_current_oi_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,30) = strpad('tau_oi in component transient_outward_K_current_oi_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,37) = strpad('oi_infinity in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,43) = strpad('g_Kur_d in component ultrarapid_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_STATES(:,9) = strpad('ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_STATES(:,10) = strpad('ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,7) = strpad('alpha_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,21) = strpad('beta_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,31) = strpad('tau_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,38) = strpad('ua_infinity in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('alpha_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,22) = strpad('beta_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,32) = strpad('tau_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,39) = strpad('ui_infinity in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_CONSTANTS(:,11) = strpad('g_Kr in component rapid_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_STATES(:,11) = strpad('xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('alpha_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,23) = strpad('beta_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,33) = strpad('tau_xr in component rapid_delayed_rectifier_K_current_xr_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,40) = strpad('xr_infinity in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_CONSTANTS(:,12) = strpad('g_Ks in component slow_delayed_rectifier_K_current (nanoS_per_picoF)');
    LEGEND_STATES(:,12) = strpad('xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,10) = strpad('alpha_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,24) = strpad('beta_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,34) = strpad('tau_xs in component slow_delayed_rectifier_K_current_xs_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,41) = strpad('xs_infinity in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('g_Ca in component sarcolemmal_Ca_current (nanoS_per_picoF)');
    LEGEND_STATES(:,13) = strpad('Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,14) = strpad('d in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_STATES(:,15) = strpad('f in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_STATES(:,16) = strpad('f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,11) = strpad('d_infinity in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('tau_d in component sarcolemmal_Ca_current_d_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,12) = strpad('f_infinity in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,26) = strpad('tau_f in component sarcolemmal_Ca_current_f_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,13) = strpad('f_Ca_infinity in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('tau_f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (millisecond)');
    LEGEND_CONSTANTS(:,15) = strpad('g_Cl_Ca in component Ca_activated_Cl_current (nanoS_per_picoF)');
    LEGEND_ALGEBRAIC(:,48) = strpad('E_Cl in component Ca_activated_Cl_current (millivolt)');
    LEGEND_ALGEBRAIC(:,58) = strpad('Fn in component Ca_release_current_from_JSR (dimensionless)');
    LEGEND_STATES(:,17) = strpad('Cl_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,16) = strpad('Cl_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,17) = strpad('q_Ca in component Ca_activated_Cl_current_q_Ca_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,60) = strpad('q_Ca_infinity in component Ca_activated_Cl_current_q_Ca_gate (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('Km_Na_i in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,19) = strpad('Km_K_o in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,20) = strpad('i_NaK_max in component sodium_potassium_pump (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,50) = strpad('f_NaK in component sodium_potassium_pump (dimensionless)');
    LEGEND_CONSTANTS(:,41) = strpad('sigma in component sodium_potassium_pump (dimensionless)');
    LEGEND_CONSTANTS(:,21) = strpad('I_NaCa_max in component Na_Ca_exchanger_current (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,22) = strpad('K_mNa in component Na_Ca_exchanger_current (millimolar)');
    LEGEND_CONSTANTS(:,23) = strpad('K_mCa in component Na_Ca_exchanger_current (millimolar)');
    LEGEND_CONSTANTS(:,24) = strpad('K_sat in component Na_Ca_exchanger_current (dimensionless)');
    LEGEND_CONSTANTS(:,25) = strpad('Ca_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,26) = strpad('g_B_Na in component background_currents (nanoS_per_picoF)');
    LEGEND_CONSTANTS(:,27) = strpad('g_B_Ca in component background_currents (nanoS_per_picoF)');
    LEGEND_ALGEBRAIC(:,54) = strpad('E_Ca in component background_currents (millivolt)');
    LEGEND_CONSTANTS(:,28) = strpad('i_p_Ca_max in component Ca_pump_current (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,57) = strpad('i_rel in component Ca_release_current_from_JSR (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,29) = strpad('K_rel in component Ca_release_current_from_JSR (per_millisecond)');
    LEGEND_CONSTANTS(:,30) = strpad('V_rel in component Ca_release_current_from_JSR (micrometre_3)');
    LEGEND_STATES(:,18) = strpad('Ca_rel in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,19) = strpad('u in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_STATES(:,20) = strpad('v in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_STATES(:,21) = strpad('w in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_CONSTANTS(:,42) = strpad('tau_u in component Ca_release_current_from_JSR_u_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,61) = strpad('u_infinity in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,62) = strpad('tau_v in component Ca_release_current_from_JSR_v_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,64) = strpad('v_infinity in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('tau_w in component Ca_release_current_from_JSR_w_gate (millisecond)');
    LEGEND_ALGEBRAIC(:,27) = strpad('w_infinity in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,59) = strpad('i_tr in component transfer_current_from_NSR_to_JSR (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,31) = strpad('tau_tr in component transfer_current_from_NSR_to_JSR (millisecond)');
    LEGEND_STATES(:,22) = strpad('Ca_up in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,32) = strpad('I_up_max in component Ca_uptake_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_ALGEBRAIC(:,63) = strpad('i_up in component Ca_uptake_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,33) = strpad('K_up in component Ca_uptake_current_by_the_NSR (millimolar)');
    LEGEND_ALGEBRAIC(:,66) = strpad('i_up_leak in component Ca_leak_current_by_the_NSR (picoA_per_picoF)');
    LEGEND_CONSTANTS(:,34) = strpad('Ca_up_max in component Ca_leak_current_by_the_NSR (millimolar)');
    LEGEND_CONSTANTS(:,35) = strpad('CMDN_max in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS(:,36) = strpad('TRPN_max in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS(:,37) = strpad('CSQN_max in component Ca_buffers (millimolar)');
    LEGEND_ALGEBRAIC(:,67) = strpad('J_Ca_CMDN in component Ca_buffers (millimolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,68) = strpad('J_Ca_TRPN in component Ca_buffers (millimolar_per_millisecond)');
    LEGEND_ALGEBRAIC(:,65) = strpad('J_Ca_CSQN in component Ca_buffers (millimolar_per_millisecond)');
    LEGEND_STATES(:,23) = strpad('Ca_CMDN in component Ca_buffers (millimolar)');
    LEGEND_STATES(:,24) = strpad('Ca_TRPN in component Ca_buffers (millimolar)');
    LEGEND_STATES(:,25) = strpad('Ca_CSQN in component Ca_buffers (millimolar)');
    LEGEND_CONSTANTS(:,38) = strpad('V_i in component intracellular_ion_concentrations (micrometre_3)');
    LEGEND_CONSTANTS(:,39) = strpad('V_rel in component intracellular_ion_concentrations (micrometre_3)');
    LEGEND_CONSTANTS(:,40) = strpad('V_up in component intracellular_ion_concentrations (micrometre_3)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,3) = strpad('d/dt m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt j in component fast_sodium_current_j_gate (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt oa in component transient_outward_K_current_oa_gate (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt oi in component transient_outward_K_current_oi_gate (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)');
    LEGEND_RATES(:,11) = strpad('d/dt xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless)');
    LEGEND_RATES(:,12) = strpad('d/dt xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless)');
    LEGEND_RATES(:,14) = strpad('d/dt d in component sarcolemmal_Ca_current_d_gate (dimensionless)');
    LEGEND_RATES(:,15) = strpad('d/dt f in component sarcolemmal_Ca_current_f_gate (dimensionless)');
    LEGEND_RATES(:,16) = strpad('d/dt f_Ca in component sarcolemmal_Ca_current_f_Ca_gate (dimensionless)');
    LEGEND_RATES(:,19) = strpad('d/dt u in component Ca_release_current_from_JSR_u_gate (dimensionless)');
    LEGEND_RATES(:,20) = strpad('d/dt v in component Ca_release_current_from_JSR_v_gate (dimensionless)');
    LEGEND_RATES(:,21) = strpad('d/dt w in component Ca_release_current_from_JSR_w_gate (dimensionless)');
    LEGEND_RATES(:,23) = strpad('d/dt Ca_CMDN in component Ca_buffers (millimolar)');
    LEGEND_RATES(:,24) = strpad('d/dt Ca_TRPN in component Ca_buffers (millimolar)');
    LEGEND_RATES(:,25) = strpad('d/dt Ca_CSQN in component Ca_buffers (millimolar)');
    LEGEND_RATES(:,2) = strpad('d/dt Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,6) = strpad('d/dt K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,17) = strpad('d/dt Cl_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,13) = strpad('d/dt Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,22) = strpad('d/dt Ca_up in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,18) = strpad('d/dt Ca_rel in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [t, STATES, ALGEBRAIC, CONSTANTS] = Kv(X, holding_p, holding_t, P1, P1_t, P2, P2_t)
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel(X, holding_p, holding_t, P1, P1_t, P2, P2_t);
end

function [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel(X, holding_p, holding_t, P1, P1_t, P2, P2_t)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts();

    % Set timespan to solve over 
    tspan = [0, P2_t];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [t, STATES] = ode15s(@(t, STATES)computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2, X), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2, X);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t, holding_p, holding_t, P1, P1_t, P2, X);
end

function [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2, X)
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
    ALGEBRAIC(:,16) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);
    
    % Ito; 10 control variables
    % A71; alpha_a
    ALGEBRAIC(:,1) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,16)+X(1)));
    % A72; beta_a
    ALGEBRAIC(:,5) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,16)+X(2)));
    % A73; alpha_i
    ALGEBRAIC(:,2) = ( 0.000152000.*exp( - (ALGEBRAIC(:,16)+X(3))./X(4)))./( 0.00670830.*exp( - (ALGEBRAIC(:,16)+X(5))./X(6))+1.00000);
    % A74; beta_i
    ALGEBRAIC(:,6) = ( 0.000950000.*exp((ALGEBRAIC(:,16)+X(7))./X(8)))./( 0.0513350.*exp((ALGEBRAIC(:,16)+X(9))./X(10))+1.00000);
    % A69; ato_f
    RATES(:,1) =  ALGEBRAIC(:,1).*(1.00000 - STATES(:,1)) -  ALGEBRAIC(:,5).*STATES(:,1);
    % A70; ito_f
    RATES(:,2) =  ALGEBRAIC(:,2).*(1.00000 - STATES(:,2)) -  ALGEBRAIC(:,6).*STATES(:,2);
    % A67; I_Kto,f
    ALGEBRAIC(:,13) =  CONSTANTS(:,1).*power(STATES(:,1), 3.00000).*STATES(:,2).*(ALGEBRAIC(:,16) + 82.8);

    % IKslow1; 6 control variables
    % A78; a_ss
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,16)+X(11))./X(12)));
    % A79; i_ss
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((ALGEBRAIC(:,16)+X(13))./X(14)));
    % A90; tau_aur
    ALGEBRAIC(:,7) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,16))+2.05800;
    % A91; tau_iur
    ALGEBRAIC(:,8) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,16)+X(15))./X(16)));
    % A88; aur
    RATES(:,3) = (ALGEBRAIC(:,3) - STATES(:,3))./ALGEBRAIC(:,7);
    % A89; iur
    RATES(:,4) = (ALGEBRAIC(:,4) - STATES(:,4))./ALGEBRAIC(:,8);
    % A87; I_kUR
    ALGEBRAIC(:,14) =  CONSTANTS(:,2).*STATES(:,3).*STATES(:,4).*(ALGEBRAIC(:,16) + 82.8);
    
    % IKslow2; 6 control variables
    % A78; a_ss
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,16)+X(17))./X(18)));
    % A79; i_ss
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((ALGEBRAIC(:,16)+X(19))./X(20)));
    % A90; tau_aur
    ALGEBRAIC(:,11) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,16))+2.05800;
    % A91; tau_iur
    ALGEBRAIC(:,12) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,16)+X(21))./X(22)));
    % A88; aur
    RATES(:,5) = (ALGEBRAIC(:,9) - STATES(:,5))./ALGEBRAIC(:,11);
    % A89; iur
    RATES(:,6) = (ALGEBRAIC(:,10) - STATES(:,6))./ALGEBRAIC(:,12);
    % A87; I_kUR
    ALGEBRAIC(:,15) =  CONSTANTS(:,2).*STATES(:,5).*STATES(:,6).*(ALGEBRAIC(:,16) + 82.8);

    RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t, holding_p, holding_t, P1, P1_t, P2, X)
    ALGEBRAIC(:,16) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);

    ALGEBRAIC(:,1) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,16)+X(1)));
    ALGEBRAIC(:,5) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,16)+X(2)));
    ALGEBRAIC(:,2) = ( 0.000152000.*exp( - (ALGEBRAIC(:,16)+X(3))./X(4)))./( 0.00670830.*exp( - (ALGEBRAIC(:,16)+X(5))./X(6))+1.00000);
    ALGEBRAIC(:,6) = ( 0.000950000.*exp((ALGEBRAIC(:,16)+X(7))./X(8)))./( 0.0513350.*exp((ALGEBRAIC(:,16)+X(9))./X(10))+1.00000);
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,16)+X(11))./X(12)));
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((ALGEBRAIC(:,16)+X(13))./X(14)));
    ALGEBRAIC(:,7) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,16))+2.05800;
    ALGEBRAIC(:,8) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,16)+X(15))./X(16)));
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,16)+X(17))./X(18)));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((ALGEBRAIC(:,16)+X(19))./X(20)));
    ALGEBRAIC(:,11) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,16))+2.05800;
    ALGEBRAIC(:,12) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,16)+X(21))./X(22)));
    % Ito
    ALGEBRAIC(:,13) =  CONSTANTS(:,1).*power(STATES(:,1), 3.00000).*STATES(:,2).*(ALGEBRAIC(:,16) + 82.8);
    % IKslow1
    ALGEBRAIC(:,14) =  CONSTANTS(:,2).*STATES(:,3).*STATES(:,4).*(ALGEBRAIC(:,16) + 82.8);
    % IKslow2
    ALGEBRAIC(:,15) =  CONSTANTS(:,2).*STATES(:,5).*STATES(:,6).*(ALGEBRAIC(:,16) + 82.8);
end

function VC = volt_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    if t < holding_t
        VC = holding_p;
    elseif (t >= holding_t) && (t <= P1_t) 
        VC = P1;
    else
        VC = P2;
    end
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % There are a total of 41 entries in each of the rate and state variable arrays.
    % There are a total of 73 entries in the constant variable array.
    algebraicVariableCount =16;
end

function [STATES, CONSTANTS] = initConsts()
    CONSTANTS = []; STATES = [];

    STATES(:,1) = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    STATES(:,2) = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    STATES(:,3) = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,4) = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    STATES(:,5) = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,6) = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    CONSTANTS(:,1) = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    CONSTANTS(:,2) = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF

    if (isempty(STATES)), warning('Initial values for states not set'); end
end

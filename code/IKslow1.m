
function [t, STATES, ALGEBRAIC, CONSTANTS] = IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, X)
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
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);
    
    % IKslow1; 6 control variables
    % A78; a_ss
    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,6)+X(1))./X(2)));
    % A79; i_ss
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((ALGEBRAIC(:,6)+X(3))./X(4)));
    % A90; tau_aur
    ALGEBRAIC(:,3) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,6))+2.05800;
    % A91; tau_iur
    ALGEBRAIC(:,4) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,6)+X(5))./X(6)));
    % A88; aur
    RATES(:,1) = (ALGEBRAIC(:,1) - STATES(:,1))./ALGEBRAIC(:,3);
    % A89; iur
    RATES(:,2) = (ALGEBRAIC(:,2) - STATES(:,2))./ALGEBRAIC(:,4);
    % A87; I_kUR
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*STATES(:,1).*STATES(:,2).*(ALGEBRAIC(:,6) + 82.8);
    
    RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t, holding_p, holding_t, P1, P1_t, P2, X)
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);

    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,6)+X(1))./X(2)));
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((ALGEBRAIC(:,6)+X(3))./X(4)));
    ALGEBRAIC(:,3) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,6))+2.05800;
    ALGEBRAIC(:,4) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,6)+X(5))./X(6)));
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*STATES(:,1).*STATES(:,2).*(ALGEBRAIC(:,6) + 82.8);
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
    algebraicVariableCount =6;
end

function [STATES, CONSTANTS] = initConsts()
    CONSTANTS = []; STATES = [];

    STATES(:,1) = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    STATES(:,2) = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    CONSTANTS(:,1) = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF

    if (isempty(STATES)), warning('Initial values for states not set'); end
end

function [t, STATES, ALGEBRAIC] = IKslow(X, hp, ht, P1, P1t, Ek)
   [t, STATES, ALGEBRAIC] = solveModel(X, hp, ht, P1, P1t, Ek);
end

function [t, STATES, ALGEBRAIC] = solveModel(X, hp, ht, P1, P1t, Ek)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  
    algebraicVariableCount = getAlgebraicVariableCount();
    
    % Initialise constants and state variables
    INIT_STATES = initConsts();

    % Set timespan to solve over 
    tspan = [0, P1t];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [t, STATES] = ode15s(@(t, STATES)computeRates(X, t, STATES, hp, ht, P1, P1t, Ek), tspan, INIT_STATES, options);
    
    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(X, t, STATES, hp, ht, P1, P1t, Ek);
    ALGEBRAIC = computeAlgebraic(X, ALGEBRAIC, STATES, t, hp, ht, P1, P1t, Ek);
end

function [RATES, ALGEBRAIC] = computeRates(X, t, STATES, hp, ht, P1, P1t, Ek)
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
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, hp, ht, P1, P1t), t);
    
    % IKslow
    % A78; a_ss
    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,6)+X(1))./X(2)));
    % A79; i_ss
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((ALGEBRAIC(:,6)+X(3))./X(4)));
    % A90; tau_aur
    ALGEBRAIC(:,3) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,6)) + 2.058;
    % A91; tau_iur
    ALGEBRAIC(:,4) = X(5) - 170.000./(1.00000+exp((ALGEBRAIC(:,6)+45.2)./5.7));
    % A88; aur
    RATES(:,1) = (ALGEBRAIC(:,1) - STATES(:,1))./ALGEBRAIC(:,3);
    % A89; iur
    RATES(:,2) = (ALGEBRAIC(:,2) - STATES(:,2))./ALGEBRAIC(:,4);
    % A87; I_kUR
    ALGEBRAIC(:,5) =  X(6).*STATES(:,1).*STATES(:,2).*(ALGEBRAIC(:,6) - Ek);
    
    RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(X, ALGEBRAIC, STATES, t, hp, ht, P1, P1t, Ek)
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, hp, ht, P1, P1t), t);

    % A78; a_ss
    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,6)+X(1))./X(2)));
    % A79; i_ss
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((ALGEBRAIC(:,6)+X(3))./X(4)));
    % A90; tau_aur
    ALGEBRAIC(:,3) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,6)) + 2.058;
    % A91; tau_iur
    ALGEBRAIC(:,4) = X(5) - 170.000./(1.00000+exp((ALGEBRAIC(:,6)+45.2)./5.7));
    % A87; I_kUR
    ALGEBRAIC(:,5) =  X(6).*STATES(:,1).*STATES(:,2).*(ALGEBRAIC(:,6) - Ek);
end

function VC = volt_clamp(t, hp, ht, P1, P1t)
    if t < ht
        VC = hp;
    elseif (t >= ht) && (t <= P1t) 
        VC = P1;
    end
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    algebraicVariableCount =6;
end

function [STATES] = initConsts()
    STATES = [];

    STATES(:,1) = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,2) = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current0
    % CONSTANTS(:,1) = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF

    if (isempty(STATES)), warning('Initial values for states not set'); end
end

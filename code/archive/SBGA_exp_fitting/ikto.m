function [current_trc] = ikto(x, hold_volt, volt, time_space, Ek)
    % default parameters
    p = [30.0, 20.0, 33.5, 7.0, 0.18064, 0.03577, 0.3956, 0.06237, 0.000152, 0.067083, 0.00095, 0.051335, 0.4067];

    % input x for selected variables
    p(1) = x(1);
    p(2) = x(2);
    p(3) = x(3);
    p(4) = x(4);
    p(13) = x(5);

    % initial values of state variables
    act0 = 0.265563e-2; % ato_f; Gating variable for transient outward K+ current
    inact0 = 0.999977; % ito_f; Gating variable for transient outward K+ current
    
    % max conductance
    gmax = p(13);

    % time space information
    t = time_space{1};
    hold_t = time_space{2};
    pulse_t = time_space{3};
    hold_idx = length(hold_t);

    current_trc = zeros(length(t), 1);

    % at holding potential
    gv_hold = gating_variables(p, hold_volt);
    act_hold = hh_model(hold_t, act0, gv_hold(1), gv_hold(2));
    inact_hold = hh_model(hold_t, inact0, gv_hold(3), gv_hold(4));
    current_trc(1:hold_idx) = gmax.*(act_hold.^3).*(inact_hold).*(hold_volt - Ek);

    % at pulse voltage
    gv_pulse = gating_variables(p, volt);
    act_pulse = hh_model(pulse_t, act0, gv_pulse(1), gv_pulse(2));
    inact_pulse = hh_model(pulse_t, inact0, gv_pulse(3), gv_pulse(4));
    current_trc((hold_idx + 1):end) = gmax.*(act_pulse.^3).*(inact_pulse).*(volt - Ek);
end

function [gv] = gating_variables(p, V)
    % gv(1): steady-state activation
    % gv(2): time constant of activation
    % gv(3): steady-state inactivation
    % gv(4): time constant of inactivation
    
    gv = zeros(4, 1);
    alphaA = p(5).*exp(p(6).*(V+p(1)));
    betaA = p(7).*exp(-p(8).*(V+p(1)));
    alphaI = (p(9).*exp(-(V+p(3)-p(2))./p(4))) ./ (p(10).*exp(-(V+p(3))./p(4)) + 1);
    betaI = (p(11).*exp((V+p(3))./p(4))) ./ (p(12).*exp((V+p(3))./p(4)) + 1);

    gv(1) = alphaA/(alphaA+betaA);
    gv(2) = 1/(alphaA+betaA);
    gv(3) = alphaI/(alphaI+betaI);
    gv(4) = 1/(alphaI+betaI);
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end

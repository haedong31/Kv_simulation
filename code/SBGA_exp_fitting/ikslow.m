function [current_trc] = ikslow(x, hold_volt, volt, time_space, Ek)
    % default parameters
    p = [22.5, 7.7, 45.2, 5.7, 0.493, 0.0629, 2.058, 1200, 170, 0.16];
    
    % input x for selected variables
    p(1) = x(1);
    p(2) = x(2);
    p(3) = x(3);
    p(4) = x(4);
    p(8) = x(5);
    p(10) = x(6);

    % initial values of state variables
    act0 = 0.417069e-3;  % aur; Gating variable for 
    inact0 = 0.998543;  % iur; Gating variable for ultrarapidly 
    
    % max conductance
    gmax = p(10); % 0.16

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
    gv = zeros(4, 1);

    % activation
    gv(1) = 1 ./ (1 + exp(-(V+p(1))./p(2)));
    gv(2) = p(5).*exp(-p(6).*V) + p(7);

    % inactivation
    gv(3) = 1 ./ (1 + exp((V+p(3))./p(4)));
    gv(4) = p(8) - (p(9))./(1 + exp((V+p(3))./p(4)));
end

function [y] = hh_model(t, ss0, ss, tau)
    y = ss - (ss - ss0).*exp(-t./tau);
end

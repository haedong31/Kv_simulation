% SSA = table2array(na_gv_ko(:,7:20));
% SSI = table2array(na_ssi_ko(:,5:20));

function e = INa_fitness(X, SSA, SSI, Gmax)
    %% SSA
    % SSA protocol
    holding_p = -100; % mV
    holding_t = 10; % ms
    P1 = -85:5:-20; % mV
    P1_t = 120 + holding_t; % ms
    P2 = -100; % mV
    P2_t = P1_t; % ms
    
    % run simulation
    As = cell(1, length(P1));
    for i=1:length(P1)
        [~, ~, A, ~] = INa(holding_p,holding_t,P1(i),P1_t,P2,P2_t,X);
        As{i} = A;
    end

    % extract peaks and Ernest potentials at peak
    peaks = zeros(1, length(P1));
    ENas = zeros(1, length(P1));
    for i=1:length(P1)
        [peak, peak_idx] = min(As{i}(:,58));
        peaks(i) = peak;
        ENas(i) = As{i}(peak_idx,57);
    end
    SSA_hat = peaks./((P1 - ENas)*Gmax);


    %% SSI
    holding_p = -100; % mV
    holding_t = 10; % ms
    P1 = -140:5:-65; % mV
    P1_t = 500 + holding_t; % ms
    P2 = -20; % mV
    P2_t = P1_t + 10; % ms

    % run simulation
    As = cell(1, length(P1));
    for i=1:length(P1)
        [~, ~, A, ~] = INa(holding_p,holding_t,P1(i),P1_t,P2,P2_t,X);
        As{i} = A;
    end
    
    % extract peaks
    peaks = zeros(1, length(P1));
    for i=1:length(P1)
        peaks(i) = min(As{i}(:,58));
    end
    Imax = min(peaks);
    SSI_hat = peaks/Imax;
    
    
    %% calculate error
    SSA_size = size(SSA);
    SSI_size = size(SSI);
    
    SSA_err = zeros(1, SSA_size(1));
    for i = 1:SSA_size(1)
        SSA_err(i) = sum((SSA_hat - SSA(i,:)).^2);
    end
    
    SSI_err = zeros(1, SSI_size(1));
    for i = 1:SSI_size(1)
        SSI_err(i) = sum((SSI_hat - SSI(i,:)).^2);
    end
    
    e = sum(SSA_err) + sum(SSI_err);
end

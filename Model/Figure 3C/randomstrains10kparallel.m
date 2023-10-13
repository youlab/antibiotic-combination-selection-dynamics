% Random strains modified to provide sensitive and resistant cell densities for 
% lowest resistant fraction along the isobole
% Calls nondimODEs.m 

close all
clear
%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 16)


nr0 = 0.2;   % initial resistant cells
ns0 = 0.2;   % initial sensitive cells
s0 = 4;    % nutrients
b0 = 0;    % Bla
antibiotics = logspace(log10(1), log10(100), 25);   % antibiotics- 25 points between 1 and 100
inhibitors = logspace(log10(0.1), log10(10), 25);   % inhibitors- 25 points between 0.1 and 10
nstrains = 10000;   % number of strains

tfinal = 100;   % final time

%% solve ODEs for single dose

t_min = 0;
t_max = 100;
time = linspace(t_min, t_max, (t_max-t_min)*6 + 1);   % time range

% variables to be updated in calculations
data = nan(nstrains, 13); % for 11 parameters + sensitive, resistant cell densities
%%
tic % start timer
parfor k = 1:nstrains   % for each strain
    
    tempdata = nan(1, 13);
    
    time = linspace(0, 100, 100*6+1)
    
    nr = zeros(length(time), length(antibiotics), length(inhibitors));
    ns = zeros(length(time), length(antibiotics), length(inhibitors));
    
    xi = rand();   % nutrient release
    kappab = rand();   % antibiotic degradation by extracellular Bla
    phimax = 5*rand();   % antibiotic degradation by resistant cells
    da = 0.02;   % natural decay of antibiotic
    ha = 4*rand() + 1;   % hill coefficient of antibiotic (1-5)
    gamma = 1.1 + 0.3*rand();   % lysis coefficient (1.1-1.4)
    db = 1+9*rand();   % degradation/inhibition of extracellular Bla (1-10)
    hi = 4*rand() + 1;   % hill coefficient of inhibitor (1-5)
    alpha = 0.75 + 0.25*rand();   % growth burden of resistant cells
    betamin = rand();    % private benefit
    c = rand();   % intracellular effectiveness

    % add randomly generated values to params
    tempdata(1) = xi;
    tempdata(2) = kappab;
    tempdata(3) = phimax;
    tempdata(4) = da;
    tempdata(5) = ha;
    tempdata(6) = gamma;
    tempdata(7) = db;
    tempdata(8) = hi;
    tempdata(9) = alpha;
    tempdata(10) = betamin;
    tempdata(11) = c;

    for i = 1:length(antibiotics)   % for each antibiotic concentration
        % initial values
        a0 = antibiotics(i);
        y0 = [nr0; ns0; 0; a0; s0];   

        for j = 1:length(inhibitors)   % for each inhibitor concentration
            thisinh = inhibitors(j);

            p= [xi;kappab;phimax;da;ha;gamma;db;hi;alpha;betamin;c;thisinh];   % parameter values (initial conditions)
            options = odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative', [1:2]);   % set acceptable error
            [time, y] = ode45(@nondimODEsvariable, time, y0, options, p);   % solve differential equations with specified inputs
            % update variables with values from solved differential equation
            nr(:, i, j) = real(y(:,1));  
            ns(:, i, j) = real(y(:,2));
        end
    end
    finalT_index = find(time >= tfinal, 1);   % find first index in time where time is greater than tfinal (100)
    nratFinalT = squeeze(nr(finalT_index, :, :));   % resistant cells at final time
    nsatFinalT = squeeze(ns(finalT_index, :, :));   % sensitive cells at final time
    ODatFinalT = nratFinalT + nsatFinalT;
    fractionRatFinalT = nratFinalT ./ ODatFinalT;
    
    ODthreshold = 10^0;   % OD threshold is 1
    ODsubthreshold_values = zeros(length(inhibitors), 2);   % [inhibitors, antibiotics], values at threshold
    ODsubthreshold_indices = zeros(length(inhibitors), 2);   % [inhibitors, antibiotics], indices of threshold

    for i = 1:length(inhibitors)   % for each inhibitor concentration
        a_index = find(ODatFinalT(:, i) < ODthreshold, 1, 'first');   % find first occurrence where final cell density < threshold
        if any(a_index)   % if threshold exists
            % update values at threshold
            ODsubthreshold_values(i, 1) = inhibitors(i);
            ODsubthreshold_values(i, 2) = antibiotics(a_index);
            % update indices at threshold
            ODsubthreshold_indices(i, 1) = i;
            ODsubthreshold_indices(i, 2) = a_index;
        end
    end
    
    % find where threshold actually starts
    subthreshold_start = 0;
    start = find(ODsubthreshold_indices(:, 2) < length(antibiotics), 1, 'first');   % find first occurrence where antibiotic index < max value
    if ~(isempty(start))   % if there is a threshold
        subthreshold_start = start;   % threshold starts at start
    end
    
    if subthreshold_start > 0   % if there is a threshold for the given strain
        OD_realthreshold = squeeze(ODsubthreshold_indices(subthreshold_start:end, :));   % find real threshold
        % weird corner case with squeeze if there's only one threshold point
        if subthreshold_start == length(inhibitors)
            OD_realthreshold = OD_realthreshold';
        end
        threshold_size = size(OD_realthreshold);   % size of threshold
        minInh = 0;   % inhibitor concentration for minimum resistant fraction along isobole
        minA = 0;   % antibiotic concentration for minimum resistant fraction along isobole
        minResFrac = 1;   % minimum resistant fraction along isobole
        for j = 1:threshold_size(1)   % traverse the threshold
            if (OD_realthreshold(j, 2)>0 && OD_realthreshold(j, 1)>0)
                resFrac = fractionRatFinalT(OD_realthreshold(j, 2), OD_realthreshold(j, 1));   % resistant fraction along threshold
                if resFrac < minResFrac   % if lowest resistant fraction, update values
                    minResFrac = resFrac;
                    minInh = OD_realthreshold(j, 1);
                    minA = OD_realthreshold(j, 2);
                end
            end
        end
    end
    
    % update values in data if a threshold exists
    if (minA~=0 && minInh~=0)
        tempdata(12) = nsatFinalT(minA, minInh);   % sensitive cell density for the lowest resistant fraction along the isobole
        tempdata(13) = nratFinalT(minA, minInh);   % resistant cell density for the lowest resistant fraction along the isobole
    end
    
    data(k, :) = tempdata;
end     

toc   % end timer

save("randomstrains.mat", "data");
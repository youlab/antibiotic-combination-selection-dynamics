%% nondimensional_simulate_USEME_logscale.m 
% Calls nondimODEs.m 
close all
clear

%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 16)

% initial conditions
nr0 = 0.2;
ns0 = 0.2;
s0 = 4;
b0 = 0;    
antibiotics = logspace(log10(1), log10(100), 25); 
inhibitors = logspace(log10(1), log10(100), 25);

% parameters used to generate figures
xi=0.8; 
kappab = 0.35; %0.3; 
phimax = 1;
da=0.02; 
ha = 3;
gamma = 1.38;     
db = 1; 
hi=2;
alpha = 0.95;
betamin = 0.9; %0.8; 
c = 0.7; %0.8;

defaultparams = [xi;kappab;phimax;da;ha;gamma;db;hi;alpha;betamin;c;inhibitors(1)];
paramnames = ["\xi", "\kappa_b", "\phi_{max}", "d_a", "h_a", "\gamma", "d_b", "h_i", "\alpha", "\beta_{min}", "c"];
    
nparams = 11;
paramranges = zeros(nparams, 2); %parameters, boundaries
paramranges(1, :) = [0 1]; 
paramranges(2, :) = [0 1];
paramranges(3, :) = [0 5];
paramranges(4, :) = [0.02 0.02];
paramranges(5, :) = [1 5];
paramranges(6, :) = [1.1 1.4];
paramranges(7, :) = [1 10];
paramranges(8, :) = [1 5];
paramranges(9, :) = [0.75 1];
paramranges(10, :) = [0.5 1];
paramranges(11, :) = [0 1];

ntestpoints = 6;
paramvals = zeros(nparams, ntestpoints); %parameters, values
for i = 1:length(paramranges)
    paramvals(i, :) = linspace(paramranges(i, 1), paramranges(i, 2), 6);
end

%% solve ODEs for single dose

t_min = 0;
t_max = 100;
time = linspace(t_min, t_max, (t_max-t_min)*6 + 1); %"every 10 minutes"

nr = zeros(nparams, ntestpoints, length(antibiotics), length(inhibitors), length(time));
ns = zeros(nparams, ntestpoints, length(antibiotics), length(inhibitors), length(time));
a = zeros(nparams, ntestpoints, length(antibiotics), length(inhibitors), length(time));
b = zeros(nparams, ntestpoints, length(antibiotics), length(inhibitors), length(time));
s = zeros(nparams, ntestpoints, length(antibiotics), length(inhibitors), length(time));


for i = 1:nparams
    tic
    for j = 1:ntestpoints
        for k = 1:length(antibiotics)
            a0 = antibiotics(k);
            y0 = [nr0; ns0; 0; a0; s0];

            for m = 1:length(inhibitors)
                p= defaultparams;
                p(i) = paramvals(i, j);
                p(12) = inhibitors(m);
              
                options = odeset('RelTol',1e-10,'AbsTol',1e-10); % reduce step size, reduce artifacts
                [time, y] = ode45(@nondimODEs2, time, y0, options, p); 
                nr(i, j, k, m, :) = y(:,1); % save resistant cell density for these concentrations
                ns(i, j, k, m, :) = y(:,2); % save sensitive cell density for these concentrations
                b(i, j, k, m, :) = y(:,3);
                a(i, j, k, m, :) = y(:,4);
                s(i, j, k, m, :) = y(:,5);
            end     
        end    
    end
    toc
end

%% Get final measurements

nr = real(nr);
ns = real(ns);
n = real(nr + ns);

% final timepoint
tfinal = 100;
% find index of tfinal
finalT_index = find(time >= tfinal, 1);

% get values at final T
nratFinalT = real(squeeze(nr(:, :, :, :, finalT_index)));
nsatFinalT = real(squeeze(ns(:, :, :, :, finalT_index)));
ODatFinalT = real(nratFinalT + nsatFinalT);
fractionRatFinalT = nratFinalT ./ ODatFinalT;


%% isoboles

ODsubthreshold_values = zeros(nparams, ntestpoints, length(inhibitors), 2); % [inhibitors antibiotics]
ODsubthreshold_indices = length(antibiotics)*ones(nparams, ntestpoints, length(inhibitors), 2);

% inhibitory 
for i = 1:nparams
    for j = 1:ntestpoints
        for k = 1:length(inhibitors)
            
            a_index = find(ODatFinalT(i, j, :, k) < 0.01, 1, 'first');
            if any(a_index)
                    ODsubthreshold_values(i, j, k, 1) = inhibitors(k);
                    ODsubthreshold_values(i, j, k, 2) = antibiotics(a_index);
                    ODsubthreshold_indices(i, j, k, 1) = k;
                    ODsubthreshold_indices(i, j, k, 2) = a_index;
            end
        end
    end
end

% figure out where real threshold starts
subthreshold_start = zeros(nparams, ntestpoints);
for i = 1:nparams
    for j = 1:ntestpoints
        start = find(ODsubthreshold_indices(i, j, :, 2) < length(antibiotics), 1, 'first');
        if ~(isempty(start))
            subthreshold_start(i, j) = start;
        end
    end
end

% get fractions along the threshold (NaN otherwise)
fraction_subthreshold = NaN(nparams, ntestpoints, length(inhibitors));

for i = 1:nparams
    for j = 1:ntestpoints
        if subthreshold_start(i, j) > 0
            OD_realthreshold = squeeze(ODsubthreshold_indices(i, j, subthreshold_start(i, j):end, :));
            % handling a weird corner case with squeeze if there's only one threshold point
            if subthreshold_start(i, j) == length(antibiotics)
                OD_realthreshold = OD_realthreshold';
            end
            threshold_size = size(OD_realthreshold);
            for k = 1:threshold_size(1)
                fraction_subthreshold(i, j, subthreshold_start(i, j) + k - 1) = fractionRatFinalT(i, j, OD_realthreshold(k, 2), OD_realthreshold(k, 1));
            end
        end
    end
end

% axis square
%% Figure S3

colors = crameri('grayC', ntestpoints+2);

for i = 1:nparams
    figure(100+i)
    for j = 1:ntestpoints
        plot(inhibitors, squeeze(fraction_subthreshold(i, j, :)), 'Color', colors(j+1, :), 'LineWidth', 5)
        hold on
    end
    plot([min(inhibitors), max(inhibitors)], [0.5, 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'LineStyle', ':')
    set(gca, 'fontsize', 24)
    set(gca, 'XScale', 'log')
    ax = gca;
    ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%.1f'))');
    leg = legend(strtrim(cellstr(num2str(squeeze(paramvals(i, :))', 2))'), 'Location', 'eastoutside');
    title(leg, paramnames(i))
    ylim([0 1])
    xlim([min(inhibitors) max(inhibitors)])
    xlabel('Inhibitor concentration');
    ylabel('Resistant fraction');
    axis square
end
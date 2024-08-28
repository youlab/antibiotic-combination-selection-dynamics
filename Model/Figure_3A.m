%% Figure_3A.m 
% Calls nondimODEs.m 

close all
clear

%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 16)

%initial conditions
nr0 = 0.2;
ns0 = 0.2;
s0 = 4; 
b0 = 0;    
antibiotics = logspace(log10(1), log10(100), 25); 

% parameters used to generate figures
xi=0.8; 
alphas = linspace(0.8, 0.95, 3);
betamins = linspace(0.7, 0.9, 3);
kappab=0.35;
phimax = 1;
da=0.02; 
ha = 3;
db = 1; 
hi=2;
gamma=1.38;     
inhibitors = logspace(log10(0.1), log10(10), 25);
cs = linspace(0.7, 0.9, 3);

%% solve ODEs for single dose

t_min = 0;
t_max = 100;
time = linspace(t_min, t_max, (t_max-t_min)*6 + 1); %"every 10 minutes"

nr = zeros(length(antibiotics), length(inhibitors), length(alphas), length(betamins), length(cs), length(time));
ns = zeros(length(antibiotics), length(inhibitors), length(alphas), length(betamins), length(cs), length(time));
a = zeros(length(antibiotics), length(inhibitors), length(alphas), length(betamins), length(cs), length(time));
b = zeros(length(antibiotics), length(inhibitors), length(alphas), length(betamins), length(cs), length(time));
s = zeros(length(antibiotics), length(inhibitors), length(alphas), length(betamins), length(cs), length(time));

for i = 1:length(antibiotics)
    a0 = antibiotics(i);
    y0 = [nr0; ns0; 0; a0; s0];
    tic
    for j = 1:length(inhibitors)
        thisinh = inhibitors(j);
        
        for k = 1:length(alphas)
            thisalpha = alphas(k);
            
            for l = 1:length(betamins)
                thisbetamin = betamins(l);

                for m = 1:length(cs)
                    thisc = cs(m);
                    
                    %[i j k l m]

                    p= [xi;kappab;phimax;da;ha;gamma;thisbetamin;thisalpha;db;hi;thisinh;thisc];
                    options = odeset('RelTol',1e-10,'AbsTol',1e-10); % reduce step size, reduce artifacts
                    [time, y] = ode45(@nondimODEs, time, y0, options, p); 
                    nr(i, j, k, l, m, :) = y(:,1); % save resistant cell density for these concentrations
                    ns(i, j, k, l, m, :) = y(:,2); % save sensitive cell density for these concentrations
                    b(i, j, k, l, m, :) = y(:,3);
                    a(i, j, k, l, m, :) = y(:,4);
                    s(i, j, k, l, m, :) = y(:,5);
                end
            end
        end
    end    
    toc
end    

%% Get final measurements

% final timepoint
tfinal = 100;
% find index of tfinal
finalT_index = find(time >= tfinal, 1);

% get values at final T
nratFinalT = squeeze(nr(:, :, :, :, :, finalT_index));
nsatFinalT = squeeze(ns(:, :, :, :, :, finalT_index));
ODatFinalT = nratFinalT + nsatFinalT;
fractionRatFinalT = nratFinalT ./ ODatFinalT;

%% Figure 3A

% alpha beta c

parameters = [3, 3, 1; %base
                1, 3, 1; %higher burden (lower alpha)
                3, 1, 1; %higher private benefit (lower betamin)
                3, 3, 3]; %higher intracellular inhibition (higher c)
for i = 1:4

    ialpha = parameters(i, 1);
    ibetamin = parameters(i, 2);
    ic = parameters(i, 3);

    % cell density
    figure(10*i)
    hold on
    imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(ODatFinalT(:, :, ialpha, ibetamin, ic)))
    colorbar
    caxis([0 4.5])
    colormap((cmocean('tempo')))
    ax = gca;
    ax.YDir = 'normal';
    ax.XTick = linspace(min(inhibitors), max(inhibitors), (length(inhibitors)/12+1));
    ax.YTick = linspace(min(antibiotics), max(antibiotics), (length(antibiotics)/12+1));
    ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%3g'))');
    ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%3g'))'); 
    %xlabel("Inhibitor")
    %ylabel("Antibiotic")
    %title("Cell Density")
    set(gca, 'fontsize', 30)
    axis square
    hold off

    % resistant fraction
    figure(10*i+1)
    hold on
    imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(fractionRatFinalT(:, :, ialpha, ibetamin, ic)))
    set(gca, 'color', [0 0 0]);
    colorbar
    caxis([0 1])
    colormap(cmocean('balance'))
    ax = gca;
    ax.YDir = 'normal';
    ax.XTick = linspace(min(inhibitors), max(inhibitors), (length(inhibitors)/12)+1);
    ax.YTick = linspace(min(antibiotics), max(antibiotics), (length(antibiotics)/12)+1);
    ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%3g'))');
    ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%3g'))');
    %xlabel("Inhibitor")
    %ylabel("Antibiotic")
    %title("Resistant Fraction")
    set(gca, 'fontsize', 30)
    axis square
    hold off
    
end
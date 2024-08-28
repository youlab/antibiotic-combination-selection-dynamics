%% nondimensional_simulate_USEME_logscale.m 
% Calls nondimODEs.m
% Outputsb  include: time courses, growth rate time courses, recovery time,
% resistance and resilience. 

close all
clear

%%
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', 16)

% initial coniditons
nr0 = 0.2;
ns0 = 0.2;
s0 = 4; 
b0 = 0;    
antibiotics = logspace(log10(1), log10(100), 25); 

% parameters used to generate figures
xi=0.8; 
alpha = 0.95;
betamins = linspace(0.5, 1, 6);
kappab=0.35;
phimax = 1;
da=0.02; 
ha = 3;
db = 1; 
hi=2;
gamma=1.38;     
inhibitors = logspace(log10(0.1), log10(10), 25);
cs = linspace(0.5, 1, 6);

%% solve ODEs for single dose

% time for simulations
t_min = 0;
t_max = 150;
time = linspace(t_min, t_max, (t_max-t_min)*6 + 1); %"every 10 minutes"

% save data
nr = zeros(length(antibiotics), length(inhibitors), length(betamins), length(cs), length(time));
ns = zeros(length(antibiotics), length(inhibitors), length(betamins), length(cs), length(time));
a = zeros(length(antibiotics), length(inhibitors), length(betamins), length(cs), length(time));
b = zeros(length(antibiotics), length(inhibitors), length(betamins), length(cs), length(time));
s = zeros(length(antibiotics), length(inhibitors), length(betamins), length(cs), length(time));

% iterate through different conditions
for i = 1:length(antibiotics)
    a0 = antibiotics(i);
    y0 = [nr0; ns0; 0; a0; s0];
    
    for j = 1:length(inhibitors)
        thisinh = inhibitors(j);
        
        for k = 1:length(betamins)
            thisbetamin = betamins(k);
            
            for m = 1:length(cs)
                thisc = cs(m);
                
                % set parameters and options
                p= [xi;kappab;phimax;da;ha;gamma;thisbetamin;alpha;db;hi;thisinh;thisc];
                options = odeset('RelTol',1e-10,'AbsTol',1e-10); % reduce step size, reduce artifacts
                % solve ODEs
                [time, y] = ode45(@nondimODEs, time, y0, options, p); 
                % save variable data
                nr(i, j, k, m, :) = y(:,1); 
                ns(i, j, k, m, :) = y(:,2); 
                b(i, j, k, m, :) = y(:,3);
                a(i, j, k, m, :) = y(:,4);
                s(i, j, k, m, :) = y(:,5);
            end
        end
    end     
end    
%% Get final measurements

% final timepoint
tfinal = 100;
% find index of tfinal
finalT_index = find(time >= tfinal, 1);

% get values at final T
nratFinalT = squeeze(nr(:, :, :, :, finalT_index));
nsatFinalT = squeeze(ns(:, :, :, :, finalT_index));
ODatFinalT = nratFinalT + nsatFinalT;
fractionRatFinalT = nratFinalT ./ ODatFinalT;
%% Figure 2c

%plotting parameters
iantibiotic = 5;
iinhibitor = 25;
ibetamin = 5; % 0.9
ic = 3; % 0.7

figure(1)
hold on
plot(time, squeeze(nr(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [104/255 156/255 154/255], 'LineWidth', 5);
plot(time, squeeze(ns(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [166/255 166/255 166/255], 'LineWidth', 5);
plot(time, squeeze(a(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [210/255 203/255 108/255], 'LineWidth', 5);
plot(time, squeeze(b(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [0 0 0], 'LineWidth', 5);
%leg = legend(["Resistant" "Sensitive" "Antibiotic" "Bla"], 'Location', 'eastoutside');
xlabel('Time');
xlim([0 20])
ylabel('Concentration');
set(gca, 'fontsize', 30)
ylim([0 4.5])
axis square
hold off
%% Figure 2c inset

%plotting parameters
iantibiotic = 5;
iinhibitor = 25;
ibetamin = 5; % 0.9
ic = 3; % 0.7

figure(2)
hold on
plot(time, squeeze(nr(iantibiotic, iinhibitor, ibetamin, ic, :)) ./ (squeeze(nr(iantibiotic, iinhibitor, ibetamin, ic, :)) + squeeze(ns(iantibiotic, iinhibitor, ibetamin, ic, :))), 'k', 'LineWidth', 10);
%xlabel('Time');
xlim([0 20])
ylim([0 1])
%ylabel('Resistant Fraction');
set(gca, 'fontsize', 48)
axis square
hold off

%% Figure 2d timecourses

%plotting parameters
iantibiotic = 25;
iinhibitor = 7;
ibetamin = 5; % 0.9
ic = 3; % 0.7

% high antibiotic low inhibitor  
figure(3)
hold on
plot(time, squeeze(nr(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [104/255 156/255 154/255], 'LineWidth', 10);
plot(time, squeeze(ns(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [166/255 166/255 166/255], 'LineWidth', 10);
%leg = legend(["Resistant" "Sensitive"], 'Location', 'eastoutside');
%xlabel('Time');
xlim([0 150])
ylim([0.0000001 10])
yticks(logspace(-7, 1, 3)) 
%ylabel('Cell density');
set(gca, 'fontsize', 48)
axis square
set(gca, 'YScale', 'log')

% low antibiotic high inhibitor  
iantibiotic = 11;
iinhibitor = 25;
ibetamin = 5; % 0.9
ic = 3; % 0.7

figure(4)
hold on
plot(time, squeeze(nr(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [104/255 156/255 154/255], 'LineWidth', 10);
plot(time, squeeze(ns(iantibiotic, iinhibitor, ibetamin, ic, :)), 'Color', [166/255 166/255 166/255], 'LineWidth', 10);
%leg = legend(["Resistant" "Sensitive"], 'Location', 'eastoutside');
%xlabel('Time');
xlim([0 150])
ylim([0.0000001 10])
yticks(logspace(-7, 1, 3)) 
%ylabel('Cell density');
set(gca, 'fontsize', 48)
set(gca, 'YScale', 'log')
axis square
hold off
%% Figure 2d heatmaps 

ibetamin = 5;
ic = 3;

% cell density 
figure(5)
hold on
imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(ODatFinalT(:, :, ibetamin, ic)))
caxis([0 4.5])
colorbar
colormap((cmocean('tempo')))
ax = gca;
ax.YDir = 'normal';
ax.XTick = linspace(min(inhibitors), max(inhibitors), (length(inhibitors)/12+1));
ax.YTick = linspace(min(antibiotics), max(antibiotics), (length(antibiotics)/12+1));
ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%2g'))');
ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%2g'))'); 
% xlabel("Inhibitor")
% ylabel("Antibiotic")
% title("Cell Density")
set(gca, 'fontsize', 30)
axis square
hold off

% resistant fraction
figure(6)
hold on
imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(fractionRatFinalT(:, :, ibetamin, ic)))
set(gca, 'color', [0 0 0]);
colorbar
caxis([0 1])
colormap(cmocean('balance'))
ax = gca;
ax.YDir = 'normal';
ax.XTick = linspace(min(inhibitors), max(inhibitors), (length(inhibitors)/12)+1);
ax.YTick = linspace(min(antibiotics), max(antibiotics), (length(antibiotics)/12)+1);
ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%2g'))');
ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%2g'))');
% xlabel("Inhibitor")
% ylabel("Antibiotic")
% title("Resistant Fraction")
set(gca, 'fontsize', 30)
axis square
hold off
%% Identify points along an inhibitory isobole (for supplemental figures)
ODsubthreshold_values = zeros(length(betamins), length(cs), length(inhibitors), 2); % [inhibitors antibiotics]
ODsubthreshold_indices = length(antibiotics)*ones(length(betamins), length(cs), length(inhibitors), 2);

% for each inhibitor value, find first antibiotic value high enough for
% cell density to be below 0.1
for i = 1:length(betamins)
    for j = 1:length(cs)
        for k = 1:length(inhibitors)
            
            a_index = find(ODatFinalT(:, k, i, j) < 0.1, 1, 'first');
            if any(a_index)
                    ODsubthreshold_values(i, j, k, 1) = inhibitors(k);
                    ODsubthreshold_values(i, j, k, 2) = antibiotics(a_index);
                    ODsubthreshold_indices(i, j, k, 1) = k;
                    ODsubthreshold_indices(i, j, k, 2) = a_index;
            end
        end
    end
end

% figure out where threshold starts
subthreshold_start = zeros(length(betamins), length(cs));
for i = 1:length(betamins)
    for j = 1:length(cs)
        start = find(ODsubthreshold_indices(i, j, :, 2) < length(antibiotics), 1, 'first');
        if ~(isempty(start))
            subthreshold_start(i, j) = start;
        end
    end
end

% get fractions along the threshold (NaN otherwise)
fraction_subthreshold = NaN(length(betamins), length(cs), length(inhibitors));

for i = 1:length(betamins)
    for j = 1:length(cs)
        if subthreshold_start(i, j) > 0
            OD_realthreshold = squeeze(ODsubthreshold_indices(i, j, subthreshold_start(i, j):end, :));
            % handling a weird corner case with squeeze if there's only one threshold point
            if subthreshold_start(i, j) == length(antibiotics)
                OD_realthreshold = OD_realthreshold';
            end
            threshold_size = size(OD_realthreshold);
            for k = 1:threshold_size(1)
                fraction_subthreshold(i, j, subthreshold_start(i, j) + k - 1) = fractionRatFinalT(OD_realthreshold(k, 2), OD_realthreshold(k, 1), i, j);
            end
        end
    end
end

%% Supplemental Figure S2

plottingsubthreshold = squeeze(ODsubthreshold_indices(ibetamin, ic, :, :));
plottingsubthreshold = cat(2, plottingsubthreshold, nan(length(inhibitors), 1));
for i = 1:length(inhibitors)
    if plottingsubthreshold(i, 2) ~= 25
        plottingsubthreshold(i, 3) = antibiotics(plottingsubthreshold(i, 2));
    end
end

% cell density
figure(7)
hold on
imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(ODatFinalT(:, :, ibetamin, ic)))
colorbar
colormap((cmocean('tempo')))
ax = gca;
ax.YDir = 'normal';
ax.XTick = [0.1 1 10];
ax.YTick = [1 10 100];
ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%2g'))');
ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%2g'))'); 
xlabel("Inhibitor")
ylabel("Antibiotic")
title("Cell Density")
set(gca, 'fontsize', 20)
scatter(inhibitors, plottingsubthreshold(:, 3), 'filled', 'k')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis square
hold off

% resistant fraction
plotfractionRatFinalT = fractionRatFinalT;
plotfractionRatFinalT(ODatFinalT<ODthreshold) = NaN;
fractionAlpha = ones(size(plotfractionRatFinalT));
fractionAlpha(isnan(plotfractionRatFinalT)) = 0;

figure(8)
hold on
imagesc([min(inhibitors) max(inhibitors)], [max(antibiotics) min(antibiotics)], flipud(fractionRatFinalT(:, :, ibetamin, ic)))
colorbar
caxis([0 1])
colormap(cmocean('balance'))
ax = gca;
ax.YDir = 'normal';
ax.XTick = [0.1 1 10];
ax.YTick = [1 10 100];
ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%2g'))');
ax.YTickLabel = strtrim(cellstr(num2str(antibiotics(1:12:end)', '%2g'))');
xlabel("Inhibitor")
ylabel("Antibiotic")
title("Resistant Fraction")
scatter(inhibitors, plottingsubthreshold(:, 3), 'filled', 'k')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'fontsize', 20)
axis square
hold off

% isobole by inhibitor
colors = cmocean('-gray', 6);

figure(9)
hold on
handles = gobjects(4, 1);
for i = 1:4
    handles(i) = plot(inhibitors, subthresholdsbytime(i, :), 'Color', colors(i+1, :), 'LineWidth', 2);
end
plot([min(inhibitors), max(inhibitors)], [0.5, 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'LineStyle', ':')
set(gca, 'fontsize', 20)
ax = gca;
ax.XTickLabel = strtrim(cellstr(num2str(inhibitors(1:12:end)', '%.1f'))');
set(gca, 'XScale', 'log')
ylim([0 1])
xlim([min(inhibitors) max(inhibitors)])
xlabel('Inhibitor concentration');
ylabel('Resistant fraction');
leg = legend(handles, string([0 20 100 150]), 'Location', 'eastoutside');
title(leg, "Time")
axis square

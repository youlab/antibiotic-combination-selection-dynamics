clear
close all

%%

nrepeats = 10;
nisolates = 311;
nparams = 10;

params = zeros(nrepeats, nisolates, 10);
guesses = zeros(nrepeats, nisolates, 10);

load('simulateddata14.mat'); % fullsimulated, simulated from avg, experimental
load('growthrateparams2.mat'); 

for i = 1:nrepeats
    params(i, :, :) = readmatrix(strcat("params_v14_", string(i-1), ".csv"), 'Range', 'B2:K312');
    guesses(i, :, :) = readmatrix(strcat("guesses_v14_", string(i-1), ".csv"), 'Range', 'B2:K312');
end

params(:, :, 1:3) = grparams;
guesses(:, :, 1:3) = grguesses;

load('anonisolatenames.mat');
isolatenames = combinedisolatenames;
parameternames = ["\mu_{max}" "K_s" "\theta" "L_n" "\kappa_b" "\phi_{max}" "\gamma" "\beta_{min}" "d_b" "c"];
%%
nonscipyfile = "resistance resilience beta.xlsx";

resistance = readmatrix(nonscipyfile, 'Range', 'A:A');
resilience = readmatrix(nonscipyfile, 'Range', 'B:B');
%% understand the variability
colors = crameri('oslo', nrepeats+2);

parameter_bounds = {[0, 2], [0, 0.4], [0, 5], [0, 0.8], [0, 6], [0, 6], [0.0, 6], [0, 1], [0, 6], [0, 1]};

figure(1)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    for j = 1:nrepeats
        scatter(1:nisolates, squeeze(guesses(j, :, i)), 15, colors(j+1, :), 'filled');
    end
    ylim(parameter_bounds{i})
    title(parameternames(i))
    set(gca, 'FontSize', 20)
end

figure(2)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    for j = 1:nrepeats
        scatter(1:nisolates, squeeze(params(j, :, i)), 15, colors(j+1, :), 'filled');
    end
    ylim(parameter_bounds{i})
    title(parameternames(i))
    set(gca, 'FontSize', 20)
end
%%
means = squeeze(mean(params, 1));
sem = squeeze(std(params, 1)) ./ sqrt(nrepeats);
averagesem = mean(sem, 1);

parameterorder = [1 2 3 8 10 7 4 5 6 9];

%individual with error
figure(3)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    errorbar(1:nisolates, means(:, parameterorder(i)), sem(:, parameterorder(i)), '.', 'MarkerSize', 10, 'Color', colors(4, :), 'CapSize', 0);
    title(parameternames(parameterorder(i)))
    ylim(parameter_bounds{parameterorder(i)})
    set(gca, 'FontSize', 16)
end

%% rsquared and plotting
colors = crameri('oslo', nrepeats+2);

r = corrcoef(simulatedfromavg, experimental);

figure(4)
for i = 1:nisolates
    for j = 1:3
        scatter(simulatedfromavg(i, j, :), experimental(i, j, :), 1, colors(4, :), 'filled')
        hold on
    end
end
plot([0 2], [0 2], 'k', 'LineWidth', 3)
axis square
xlim([0 2])
yticks([0 1 2])
ylim([0 2])
set(gca, 'FontSize', 24)
xlabel("Simulated")
ylabel('Experimental')
text(1.9, 0.2, strcat("r = ", compose('%.3g', r(1, 2))), 'HorizontalAlignment', 'right', 'FontSize', 30)

%% simulated vs experimental

colors = [254 224 144;
          215 48 39;
          69 117 180];
colors = colors ./ 255;

timepoints = repmat(linspace(0, 24, 145), 3, 1);
for i = 1:3
    timepoints(i, :) = 24*(i - 1) + timepoints(i, :);
end

nrows = 10;
ncols = 8;
plotsperfig = nrows * ncols;
maxfigs = ceil(nisolates / plotsperfig);
for k = 1:maxfigs
    f = figure(100+k);
    offset = plotsperfig*(k - 1);
    if k == maxfigs
        jrange = offset+1:nisolates;
    else
        jrange = (1:plotsperfig) + offset;
    end 
    t = tiledlayout(nrows, ncols);
    for j = jrange
        nexttile
        %subplot(nrows, ncols, j - offset)
        hold on
        for i = 1:3
            plot(timepoints(i, :), squeeze(simulatedfromavg(j, i, :)), 'Color', [0.7 0.7 0.7], 'LineWidth', 3)
            plot(timepoints(i, :), squeeze(experimental(j, i, :)), 'Color', colors(i, :), 'LineWidth', 3)
        end
        xlim([0 max(timepoints, [], 'all')])
        ylim([0 2])
        set(gca, 'xTick', [], 'ytick', [], 'XTickLabels', [], 'YTickLabels', []);
        text(1, 1.8, isolatenames(j), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
    end
    t.TileSpacing = 'compact';
    %t.Padding = 'compact';
    f.WindowState = 'maximized';
    %saveas(f, strcat("newsimvsexp", string(k)), 'jpeg')
end
%% mumax vs betamin
betamins = means(:, 8);
experimentalisolatenames = ["2430", "2620", "2744", "4682", "5833", "5865", "6156", "D021"];
experimentalisolates = find(matches(isolatenames, experimentalisolatenames));

figure(6)
scatter(squeeze(means(:, 1)), betamins, 50, 'k', 'filled')
hold on
scatter(squeeze(means(experimentalisolates, 1)), betamins(experimentalisolates), 50, 'r', 'filled')
axis square
set(gca, 'FontSize', 30)
ylim([0 1])
xlim(parameter_bounds{1})
xlabel('\mu_{max}')
ylabel('\beta_{min}')
%% betamin x resilience

figure(7)
scatter(resilience, betamins, 50, 'k', 'filled')
hold on
scatter(resilience(experimentalisolates), betamins(experimentalisolates), 50, 'r', 'filled')
axis square
set(gca, 'FontSize', 30)
xlabel("resilience")
ylabel("\beta_{min}")

%% betamin x resistance

figure(8)
scatter(resistance, betamins, 50, 'k', 'filled')
hold on
scatter(resistance(experimentalisolates), betamins(experimentalisolates), 50, 'r', 'filled')
axis square
set(gca, 'FontSize', 30)
xlabel("resistance")
ylabel("\beta_{min}")

%% histogram of resistant-only betamins

colors = cmocean('amp', 10);

resistantonlyisolates = find(resilience > 0);
sensitiveisolates = find(resilience <= 0);
figure(9)
hold on
histogram(betamins(resistantonlyisolates), [0:0.05:1], 'FaceColor', colors(9, :))
histogram(betamins(sensitiveisolates), [0:0.05:1], 'FaceColor', colors(2, :))
axis square
xlim([0 1])
ylim([0 100])
set(gca, 'FontSize', 26)
xlabel('estimated \beta_{min}')
ylabel('number of isolates')

%% save means
avgparams = means;
%save("estimated parameters.mat", "avgparams", "parameternames", "combinedisolatenames");
clear
close all

%%

nrepeats = 10;
nisolates = 311;
nparams = 10;

uniformparams = zeros(nrepeats, nisolates, 10);
uniformguesses = zeros(nrepeats, nisolates, 10);
gaussianparams = zeros(nrepeats, nisolates, 10);
gaussianguesses = zeros(nrepeats, nisolates, 10);

load('simulateddata.mat'); % fullsimulated, simulated from avg, experimental

for i = 1:nrepeats
    uniformparams(i, :, :) = readmatrix(strcat("params_v13_", string(i-1), ".csv"), 'Range', 'B2:K312');
    uniformguesses(i, :, :) = readmatrix(strcat("guesses_v13_", string(i-1), ".csv"), 'Range', 'B2:K312');
    gaussianparams(i, :, :) = readmatrix(strcat("growthrates_v2_", string(i-1), ".csv"), 'Range', 'B2:K312');
    gaussianguesses(i, :, :) = readmatrix(strcat("growthrateguesses_v2_", string(i-1), ".csv"), 'Range', 'B2:K312');
end

load('anonisolatenames.mat');
isolatenames = combinedisolatenames;
% isolatenames = readmatrix("params_v7_0.csv", 'Range', 'A1:A312', 'OutputType', 'string');
% isolatenames = ["None" "Low Copy Bla" "Low Copy BlaM" "High Copy Bla" "High Copy BlaM"];
parameternames = ["\mu_{max}" "K_s" "\theta" "Ln" "\kappa_b" "\phi_{max}" "\gamma" "\beta_{min}" "d_b" "c"];

%% understand the variability
colors = crameri('oslo', nrepeats+2);

figure(1)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    for j = 1:nrepeats
        scatter(1:nisolates, squeeze(gaussianparams(j, :, i)), 15, colors(j+1, :), 'filled');
    end
    title(parameternames(i))
end

figure(10)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    for j = 1:nrepeats
        scatter(1:nisolates, squeeze(uniformparams(j, :, i)), 15, colors(j+1, :), 'filled');
    end
    title(parameternames(i))
end
%%
uniformmeans = squeeze(mean(uniformparams, 1));
gaussianmeans = squeeze(mean(gaussianparams, 1));
uniformsem = squeeze(std(uniformparams, 1)) ./ sqrt(nrepeats);
gaussiansem = squeeze(std(gaussianparams, 1)) ./ sqrt(nrepeats);
uniformaveragesem = mean(uniformsem, 1);
gaussianaveragsem = mean(gaussiansem, 1);

%individual with error
figure(2)
for i = 1:nparams
    subplot(2, 5, i)
    hold on
    errorbar(1:nisolates, gaussianmeans(:, i), gaussiansem(:, i), '.', 'MarkerSize', 10, 'Color', colors(4, :), 'CapSize', 0);
    title(parameternames(i))
    set(gca, 'FontSize', 16)
end

grparams = gaussianparams(:, :, 1:3);
grguesses = gaussianguesses(:, :, 1:3);
grmeans = gaussianmeans(:, 1:3);
grsem = gaussiansem(:, 1:3);

save("growthrateparams2.mat", "grmeans", "grsem", "grparams", "grguesses");
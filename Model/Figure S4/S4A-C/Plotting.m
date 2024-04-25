clear
close all

%%
load('randomstrains_10k_s2.mat') % load simulated data

%xi, kappab, phimax, da, ha, gamma, db, hi, alpha, betmain, c, sensitive density,
%resistant density
% latter two are calculated at minimum resistant fraction along the isobole

nstrains = length(data); 

maxfraccriterion = zeros(nstrains, 3);
colors = maxfraccriterion; % colors is also this size

for i = 1:nstrains
    maxfraccriterion(i, 1) = (1 - data(i, 11)) * (1 - data(i, 10)); %(1-bmin)(1-c)
    maxfraccriterion(i, 2) = (1 - data(i, 9)) / data(i, 6); %(1-alpha)/gamma
    maxfraccriterion(i, 3) = (data(i, 13) / (data(i, 13) + data(i, 12))); %resistant fraction
end

unplottedcount=0;
map = cmocean('balance', 1000, 'negative');

% color each point by resistant fraction as long as an isobole point was
% found
for i = 1:nstrains
        if ~isnan(maxfraccriterion(i, 3)) 
        colors(i, :) = map(100 + ceil(maxfraccriterion(i, 3) * 800), :);
    else
        colors(i, :) = [1 1 1];
        unplottedcount = unplottedcount+1;
    end
end
nanindices = find(isnan(maxfraccriterion(:, 3)));

%% Figure 3C

figure(1)
hold on
scatter(maxfraccriterion(:, 2), maxfraccriterion(:, 1), 10, colors, 'filled')
scatter(maxfraccriterion(nanindices, 2), maxfraccriterion(nanindices, 1), 10, [0.8 0.8 0.8]) 
plot([0 1], [0 1], 'k', 'LineWidth', 5)
colormap(map)
clim([-0.125 1.125])
%colorbar
yticks([0 0.1 0.2])
xlim([0 0.2])
ylim([0 0.2])
%ylabel('(1 - c)(1 - \betamin)');
%xlabel('(1 - \alpha ) / \gamma');
set(gca, 'fontsize', 30)
axis square

%% alternate plotting 

figure(2)
scatter((maxfraccriterion(:, 1) ./ maxfraccriterion(:, 2)), maxfraccriterion(:, 3), 10, 'k', 'filled');
ylabel('resistant fraction');
xlabel('((1 - c)(1 - \betamin))/((1 - \alpha ) / \gamma)');
ylim([0 1])
xlim([0 2.5])
set(gca, 'fontsize', 24)
axis square
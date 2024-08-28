close all
clear

load("estimated parameters.mat");
datafile = 'Isolate DA plating summary.xlsx';

opts = detectImportOptions(datafile, 'Sheet', 1, 'Range', 'A:G'); 
opts.VariableTypes = {'double', 'string', 'string', 'double', 'double', 'double', 'double'};
allinfo = readtable(datafile, opts);
alldata = readmatrix(datafile, 'Sheet', 1, 'Range', 'H:I', 'NumHeaderLines', 1); %LB LB+CM

rinfo = string(allinfo{:, 2});
rinfo = strrep(rinfo, "DICON", "D");
allinfo(:, 2) = cellstr(rinfo);

%% subset CM-treated, AMX-treated data
cmtreatedinfo = allinfo(table2array(allinfo(:, 4)) > 0, :);
cmtreateddata = alldata(table2array(allinfo(:, 4)) > 0, :);

strainnames = string(unique(cmtreatedinfo{:, 2}));

amxtreatedinfo = cmtreatedinfo(table2array(cmtreatedinfo(:, 5)) > 0, :);
amxtreateddata = cmtreateddata(table2array(cmtreatedinfo(:, 5)) > 0, :);
inhconcentrations = unique(amxtreatedinfo{:, 6});
nreplicates = max(amxtreatedinfo{:, 7});

datasplit = zeros(length(strainnames), length(inhconcentrations), nreplicates, 2);
for i = 1:length(amxtreateddata)
    datasplit(amxtreatedinfo{i, 2} == strainnames, amxtreatedinfo{i, 6} == inhconcentrations, amxtreatedinfo{i, 7}, :) = amxtreateddata(i, :);
end

%% get fractions and averages
LB = squeeze(datasplit(:, :, :, 1));
LBcm = squeeze(datasplit(:, :, :, 2));

fractions = (LB - LBcm) ./ LB;
fractionsav = mean(fractions, 3);
fractionsSEM = std(fractions, 0, 3) ./ sqrt(nreplicates);
%% betamin processing
estimatedbetamins = avgparams(:, 8);
strainindices = find(contains(combinedisolatenames, strainnames));
estimatedbetaminssubset = estimatedbetamins(strainindices);
betaminstring = compose("%.2g", estimatedbetaminssubset);
%% plots 

% for i = 1:length(strainnames)
%     f = figure(i);
%     hold on
% 
%     %plot averages
%     errorbar(inhconcentrations, fractionsav(i, :), fractionsSEM(i, :), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
% 
%     [~, minindex] = min(estimatedbetaminssubset);
% 
%     text(3.4, 0.9, strainnames(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'FontWeight', 'bold')
%     if i == minindex
%         text(3.4, 0.78, strcat("\betamin = ", betaminstring(i)), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
%     else
%         text(3.4, 0.78, betaminstring(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
%     end 
%     %title(strainnames(i))
%     axis([0 4 0 1])
%     set(gca,'FontSize', 30)
%     %xlabel("CLA");
%     %ylabel("Resistant fraction");
%     axis square
% 
%     %saveas(f, strcat(strainnames(i), " no labels"), 'svg')
% end
%% plot replicates
for i = 1:length(strainnames)
    f = figure(i);
    hold on
    % scatter replicates
    scatter(inhconcentrations, squeeze(fractions(i, :, :)), 100, [0.7 0.7 0.7], 'filled');

    %plot averages
    errorbar(inhconcentrations, fractionsav(i, :), fractionsSEM(i, :), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
    
    [~, minindex] = min(estimatedbetaminssubset);
    
    text(3.4, 0.9, strainnames(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'FontWeight', 'bold')
    if i == minindex
        text(3.4, 0.78, strcat("\betamin = ", betaminstring(i)), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
    else
        text(3.4, 0.78, betaminstring(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
    end 
    %title(strainnames(i))
    axis([0 4 0 1])
    set(gca,'FontSize', 30)
    %xlabel("CLA");
    %ylabel("Resistant fraction");
    axis square
    box on
    
    %saveas(f, strcat("Figure 5E ", strainnames(i)), 'svg')
end
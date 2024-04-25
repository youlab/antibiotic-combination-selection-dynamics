% updated to correct for issues noted in plating_summary.m

close all
clear

load("estimated parameters.mat");
datafile = 'Isolate Keio plating summary.xlsx';
opts = detectImportOptions(datafile, 'Sheet', 1, 'Range', 'A:F'); % in the future can choose to import replicates from sheet 2
opts.VariableTypes = {'double', 'string', 'string', 'double', 'double', 'double'};
allinfo = readtable(datafile, opts);
alldata = readmatrix(datafile, 'Sheet', 1, 'Range', 'G:I', 'NumHeaderLines', 1); %LB %Fraction % FractionSEM
%%
cmtreatedinfo = allinfo(table2array(allinfo(:, 4)) > 0, :);
cmtreateddata = alldata(table2array(allinfo(:, 4)) > 0, :);
strainnames = string(unique(cmtreatedinfo{:, 2}));
%%
estimatedbetamins = avgparams(:, 8);
strainindices = find(contains(combinedisolatenames, strainnames));
estimatedbetaminssubset = estimatedbetamins(strainindices);
betaminstring = compose("%.2g", estimatedbetaminssubset);
[~, minindex] = min(estimatedbetaminssubset);

% plot fractions
for i = 1:length(strainnames)
    f = figure(i);
    indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
    thisinfo = cmtreatedinfo(indices, :);
    thisdata = cmtreateddata(indices, :);
    errorbar(0:2:4, squeeze(thisdata(2:end, 2)), squeeze(thisdata(2:end, 3)), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
        
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
    
    %saveas(f, strcat(strainnames(i), " no labels"), 'jpg')
end
%%
 % plot cfu
for i = 1:length(strainnames)
    f = figure(10 + i);
    indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
    thisinfo = cmtreatedinfo(indices, :);
    thisdata = cmtreateddata(indices, :);
    plot(0:2:4, squeeze(thisdata(2:end, 1)),  '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 

    text(3.4, 5*10^9, strainnames(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'FontWeight', 'bold')
    if i == minindex
        text(3.4, 2.2*10^9, strcat("\betamin = ", betaminstring(i)), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
    else
        text(3.4, 2.2*10^9, betaminstring(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
    end 
    %title(strainnames(i))
    xlim([0 4])
    set(gca,'FontSize',30)
    set(gca, 'YScale', 'log')
    ylim([10^7 10^10])
    % xlabel("CLA");
    % ylabel("OD600");
    axis square

    saveas(f, strcat(strainnames(i), " cfu"), 'jpg')
end



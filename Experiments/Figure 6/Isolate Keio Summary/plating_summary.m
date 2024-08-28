% updated to correct

close all
clear

load("estimated parameters.mat");
datafile = 'Isolate Keio plating summary.xlsx';
opts = detectImportOptions(datafile, 'Sheet', 1, 'Range', 'A:F'); % in the future can choose to import replicates from sheet 2
opts.VariableTypes = {'double', 'string', 'string', 'double', 'double', 'double'};
allinfo = readtable(datafile, opts);
alldata = readmatrix(datafile, 'Sheet', 1, 'Range', 'G:J', 'NumHeaderLines', 1); %LB %Fraction %propagated error %FractionSEM
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
% for i = 1:length(strainnames)
%     f = figure(i);
%     indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
%     thisinfo = cmtreatedinfo(indices, :);
%     thisdata = cmtreateddata(indices, :);
%     errorbar(0:2:4, squeeze(thisdata(2:end, 2)), squeeze(thisdata(2:end, 4)), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
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
%     %saveas(f, strcat(strainnames(i), " no labels"), 'jpg')
% end

%% plot replicates

opts = detectImportOptions(datafile, 'Sheet', 2, 'Range', 'A:G'); % 
opts.VariableTypes = {'double', 'string', 'string', 'double', 'double', 'double', 'double'};
replicatesinfo = readtable(datafile, opts);
replicatesdata = readmatrix(datafile, 'Sheet', 2, 'Range', 'H:M', 'NumHeaderLines', 1); %LB %LB error %CM %cm error %Fraction % FractionSEM

cmtreatedinforeps = replicatesinfo(table2array(replicatesinfo(:, 5)) > 0, :);
cmtreateddatareps = replicatesdata(table2array(replicatesinfo(:, 5)) > 0, :);

% plot fractions
for i = 1:length(strainnames)
    f = figure(100+i);
    hold on

    indices = string(cmtreatedinforeps{:, 2}) == strainnames(i) & (cmtreatedinforeps{:, 6} == 2);
    thisinfo = cmtreatedinforeps(indices, :);
    thisdata = cmtreateddatareps(indices, :);
    thisavindices = string(cmtreatedinfo{:, 2}) == strainnames(i);
    thisavdata = cmtreateddata(thisavindices, :);

    % scatter replicates
    scatter(thisinfo{:, 7}, thisdata(:, 5), 100, [0.7 0.7 0.7], 'filled');

    % plot averages
    errorbar(0:2:4, squeeze(thisavdata(2:end, 2)), squeeze(thisavdata(2:end, 4)), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 


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
    
    saveas(f, strcat("Figure 6 ", strainnames(i)), 'svg')
end

%%
 % plot cfu
% for i = 1:length(strainnames)
%     f = figure(10 + i);
%     indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
%     thisinfo = cmtreatedinfo(indices, :);
%     thisdata = cmtreateddata(indices, :);
%     plot(0:2:4, squeeze(thisdata(2:end, 1)),  '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
% 
%     text(3.4, 5*10^9, strainnames(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'FontWeight', 'bold')
%     if i == minindex
%         text(3.4, 2.2*10^9, strcat("\betamin = ", betaminstring(i)), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
%     else
%         text(3.4, 2.2*10^9, betaminstring(i), 'HorizontalAlignment', 'right', 'FontSize', 30, 'Color', [0.7748 0.4149 0.3004])
%     end 
%     %title(strainnames(i))
%     xlim([0 4])
%     set(gca,'FontSize',30)
%     set(gca, 'YScale', 'log')
%     ylim([10^7 10^10])
%     % xlabel("CLA");
%     % ylabel("OD600");
%     axis square
% 
%     % saveas(f, strcat(strainnames(i), " cfu"), 'jpg')
% end



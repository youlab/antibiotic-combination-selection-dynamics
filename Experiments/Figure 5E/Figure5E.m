close all
clear

load("estimated parameters.mat");
datafile = 'Isolate DA plating summary.xlsx';
allinfo = readtable(datafile, 'Sheet', 1, 'Range', 'A:F');
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'G:K'); % OD LB LBSEM LB+CM LB+CMSEM

% remove identifiers
rinfo = string(allinfo{:, 2});
rinfo = strrep(rinfo, "DICON", "D");
allinfo(:, 2) = cellstr(rinfo);
%% data
fraction = (rawdata(:, 2) - rawdata(:, 4)) ./ rawdata(:, 2);
fractionsem = sqrt((rawdata(:, 3)./rawdata(:, 2)).^2 + (sqrt((rawdata(:, 5).^2 + rawdata(:, 3).^2))./(rawdata(:, 2) - rawdata(:, 4) ./ rawdata(:, 2))).^2);
alldata = cat(2, rawdata, fraction, fractionsem);
%% data with cm
cmtreatedinfo = allinfo(table2array(allinfo(:, 4)) > 0, :);
cmtreateddata = alldata(table2array(allinfo(:, 4)) > 0, :);
strainnames = string(unique(cmtreatedinfo{:, 2}));
%% plots 
estimatedbetamins = avgparams(:, 8);
strainindices = find(contains(combinedisolatenames, strainnames));
estimatedbetaminssubset = estimatedbetamins(strainindices);
betaminstring = compose("%.2g", estimatedbetaminssubset);

for i = 1:length(strainnames)
    f = figure(i);
    indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
    thisinfo = cmtreatedinfo(indices, :);
    thisdata = cmtreateddata(indices, :);
    errorbar(0:2:4, squeeze(thisdata(2:end, 6)), squeeze(thisdata(2:end, 7)), '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
    
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
    
    %saveas(f, strcat(strainnames(i), " no labels"), 'jpg')
end

%% plot ODs
for i = 1:length(strainnames)
    f = figure(10 + i);
    indices = string(cmtreatedinfo{:, 2}) == strainnames(i);
    thisinfo = cmtreatedinfo(indices, :);
    thisdata = cmtreateddata(indices, :);
    plot(0:2:4, squeeze(thisdata(2:end, 1)),  '-ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2); 
    
    title(strainnames(i))
    axis([0 4 0 2])
    set(gca,'FontSize',20)
    xlabel("CLA");
    ylabel("OD600");
    axis square
    
    %saveas(f, strcat(strainnames(i), " OD"), 'jpg')
end
close all
clear

% set constants
nwells = 96;
height = 8;
width = 12;
ntimepoints = 145;
nstrains = 11;

strainnames = ["2683", "3599", "2529", "2744", "5865", "D021", "D036", "2430", "5833", "4682", "DA28102"];
medianames = ["CM" "AMX" "CLA"];

% read in raw data from excel
datafile = '10-9-2022 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D54:CU198')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B54:B198');

layoutfile = '10-9-2022 Plate Layout.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'G1:R8');
layoutlong = reshape(platelayout', [nwells 1]);
conditioninfo = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'A:E');

cm = unique(conditioninfo(:, 3));
amx = unique(conditioninfo(:, 4));
cla = unique(conditioninfo(:, 5));
nconditions = max(unique(conditioninfo(:, 1)));

timepoints = timedata ./ 3600;

%% process data

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% smooth
datasmoothed = movmean(rawdata, 5, 2);

% blanking smoothed data
blankdata = datasmoothed(layoutlong == 0, :);
% threshold the blanks to remove contamination
cleanblankdata = blankdata(max(blankdata, [], 2) < 0.09, :);
blanks = mean(cleanblankdata);
%blanks = 0.05*ones(1, ntimepoints);
% blanks = zeros(1, ntimepoints);
blankeddata = datasmoothed - blanks;
%blankeddata(blankeddata<0) = 0;
%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

straindata = zeros(nstrains, length(cm), length(amx), length(cla), ntimepoints);

for i = 1:nconditions
    straindata(conditioninfo(i, 2), cm == conditioninfo(i, 3), amx == conditioninfo(i, 4), cla == conditioninfo(i, 5), :) = conditionsdata{i+1}(:, :);
end
%% GRs

GRs = gradient(real(log(blankeddata)), 10/60); % real

GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

strainGRs = zeros(nstrains, length(cm), length(amx), length(cla), ntimepoints);
for i = 1:nconditions
    strainGRs(conditioninfo(i, 2), cm == conditioninfo(i, 3), amx == conditioninfo(i, 4), cla == conditioninfo(i, 5), :) = GRconditions{i+1}(:, :);
end

%% LB only
cmap = crameri('romaO', nstrains-1);

figure(1)
hold on
for i = 1:nstrains - 1
    plot(timepoints, squeeze(straindata(i, 1, 1, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(straindata(nstrains, 1, 1, 1, :)), 'k:', 'LineWidth', 5)        
xlabel("Hours");
ylabel("OD");
xlim([0 24]);
ylim([0 2]);
legend(strainnames, 'Location', 'eastoutside');
set(gca,'FontSize',20)
axis square
hold off

%GRs
figure(2)
hold on
for i = 1:nstrains - 1
    plot(timepoints, squeeze(strainGRs(i, 1, 1, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(strainGRs(nstrains, 1, 1, 1, :)), 'k:', 'LineWidth', 5)        
xlabel("Hours");
ylabel("GR (1/h)");
xlim([0 10]);
ylim([-0.5 3]);
legend(strainnames, 'Location', 'eastoutside');
set(gca,'FontSize',20)
axis square
hold off

%% LB with cm

figure(3)
hold on
for i = 1:nstrains - 1
    plot(timepoints, squeeze(straindata(i, 2, 1, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(straindata(nstrains, 2, 1, 1, :)), 'k:', 'LineWidth', 5)        
xlabel("Hours");
ylabel("OD");
xlim([0 24]);
ylim([0 2]);
legend(strainnames, 'Location', 'eastoutside');
set(gca,'FontSize',20)
axis square
hold off

%GRs
figure(4)
hold on
for i = 1:nstrains - 1
    plot(timepoints, squeeze(strainGRs(i, 2, 1, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(strainGRs(nstrains, 2, 1, 1, :)), 'k:', 'LineWidth', 5)        
xlabel("Hours");
ylabel("GR (1/h)");
xlim([0 10]);
ylim([-0.5 3]);
legend(strainnames, 'Location', 'eastoutside');
set(gca,'FontSize',20)
axis square
hold off

%% strains with abx

cmap = crameri('lajolla', length(amx)*length(cla)+2);
for i = 1:length(amx)
    for j = 1:length(cla)
        medialeg(length(cla)*(i-1)+j) = strcat(amx(i) + " \mug/mL AMX + " + cla(j) + " \mug/mL CLA");
    end
end

for i = 1:nstrains
    figure(5)
    subplot(2, ceil(nstrains /2), i)
    hold on
    for j = 1:length(amx)
        for k = 1:length(cla)
            plot(timepoints, squeeze(straindata(i, 1, j, k, :)), 'Color', cmap(length(cla)*(j-1)+k+1, :), 'LineWidth', 5)  
        end
    end
    xlabel("Hours");
    ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    title(strainnames(i))
    %legend(medialeg, 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off
end

%% strains with abx cm

for i = 1:nstrains
    figure(6)
    subplot(2, ceil(nstrains /2), i)
    hold on
    for j = 1:length(amx)
        for k = 1:length(cla)
            plot(timepoints, squeeze(straindata(i, 2, j, k, :)), 'Color', cmap(length(cla)*(j-1)+k+1, :), 'LineWidth', 5)  
        end
    end
    xlabel("Hours");
    ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    title(strainnames(i))
    %legend(medialeg, 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off
end

%% experimental strains with and without abx cm
strainorder = [8 9 10 5 6 4 11];

figure(100)
for i = 1:length(strainorder)
    subplot(1, length(strainorder), i)
    hold on
    for j = 1:length(amx)
        for k = 1:length(cla)
            plot(timepoints, squeeze(straindata(strainorder(i), 1, j, k, :)), 'Color', cmap(length(cla)*(j-1)+k+1, :), 'LineWidth', 5)  
        end
    end
    if i == 1
        xlabel("Hours");
        ylabel("OD");
        set(gca, 'xtick', 0:12:24)
    else
        set(gca, 'xtick', [], 'ytick', [], 'XTickLabels', [], 'YTickLabels', [])
    end
    xlim([0 24]);
    ylim([0 2]);
    title(strainnames(strainorder(i)))
    %legend(medialeg, 'Location', 'southoutside');
    %legend('boxoff')
    set(gca,'FontSize',20)
    axis square
    hold off
end

figure(101)
for i = 1:length(strainorder)
    subplot(1, length(strainorder), i)
    hold on
    for j = 1:length(amx)
        for k = 1:length(cla)
            plot(timepoints, squeeze(straindata(strainorder(i), 2, j, k, :)), 'Color', cmap(length(cla)*(j-1)+k+1, :), 'LineWidth', 5)  
        end
    end
    if i == 1
        xlabel("Hours");
        ylabel("OD");
        set(gca, 'xtick', 0:12:24)
    else
        set(gca, 'xtick', [], 'ytick', [], 'XTickLabels', [], 'YTickLabels', [])
    end
    xlim([0 24]);
    ylim([0 2]);
    title(strainnames(strainorder(i)))
    %legend(medialeg, 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off
end
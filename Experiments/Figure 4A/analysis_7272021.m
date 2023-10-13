close all
clear

% set constants
nwells = 96;
height = 8;
width = 12;
ntimepoints = 145;
nstrains = 6;
nmedia = 2;
nconditions = nstrains * nmedia;
nreplicates = 6;

strainnames = {'MG1655' 'WT' 'Bla' 'BlaM' 'High Bla' 'High BlaM'};
medianames = {'LB' 'LB + 100 \mug/mL carb'};

% read in raw data from excel
datafile = '7-27-2021 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D92:CU236')';
GFPdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D240:CU384')';
BFPdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D388:CU532')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B92:B236');

layoutfile = '7-27-2021 Plate Layout.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'E1:P8');
layoutlong = reshape(platelayout', [nwells 1]);
conditioninfo = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'A2:C13');

timepoints = timedata ./ 3600;

%% process data

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% smooth
datasmoothed = movmean(rawdata, 5, 2);

% blanking smoothed data
blankdata = datasmoothed(layoutlong == 0, :);
% threshold the blanks to remove contamination
cleanblankdata = blankdata(max(blankdata, [], 2) < 0.07, :);
blanks = mean(cleanblankdata);
%blanks = 0.05*ones(1, ntimepoints);
% blanks = zeros(1, ntimepoints);
blankeddata = datasmoothed - blanks;

% process GFP data
GFPdatasmooth = movmean(GFPdata, 5, 2);
GFPblankdata = GFPdatasmooth(layoutlong==0, :);
cleanGFPblankdata = GFPblankdata(blankdata(:, ntimepoints) < 0.1, :);
GFPblanks = mean(cleanGFPblankdata);
GFPblankeddata = GFPdatasmooth - GFPblanks;

% process BFP data
BFPdatasmooth = movmean(BFPdata, 5, 2);
BFPblankdata = BFPdatasmooth(layoutlong==0, :);
cleanBFPblankdata = BFPblankdata(blankdata(:, ntimepoints) < 0.1, :);
BFPblanks = mean(cleanBFPblankdata);
BFPblankeddata = BFPdatasmooth - BFPblanks;
%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

GFPconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    GFPconditionsdata{i+1} = GFPblankeddata(layoutlong==i, :);
end

BFPconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    BFPconditionsdata{i+1} = BFPblankeddata(layoutlong==i, :);
end

straindata = zeros(nstrains, nmedia, nreplicates, ntimepoints);
GFPstraindata = zeros(nstrains, nmedia, nreplicates, ntimepoints);
BFPstraindata = zeros(nstrains, nmedia, nreplicates, ntimepoints);

for i = 1:nconditions
    straindata(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = conditionsdata{i+1}(:, :);
    GFPstraindata(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = GFPconditionsdata{i+1}(:, :);
    BFPstraindata(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = BFPconditionsdata{i+1}(:, :);
end
%%
% plot all replicates
% figure(1)
% hold on
% for j = 1:nmedia
%     for i = 1:nstrains
%         subplot(nmedia, nstrains, (j-1)*nstrains+i)
%         plot(timepoints, squeeze(straindata(i, j, :, :)), 'k');
%         axis([0 24 0 1.5])
%     end
% end
% 
% figure(2)
% hold on
% for j = 1:nmedia
%     for i = 1:nstrains
%         subplot(nmedia, nstrains, (j-1)*nstrains+i)
%         plot(timepoints, squeeze(GFPstraindata(i, j, :, :)), 'k');
%         axis([0 24 0 3*10^4])
%     end
% end
% 
% 
% figure(3)
% hold on
% counter = 1;
% for j = 1:nmedia
%     for i = 1:nstrains
%         subplot(nmedia, nstrains, (j-1)*nstrains+i)
%         plot(timepoints, squeeze(BFPstraindata(i, j, :, :)), 'k');
%         counter = counter+1;
%         axis([0 24 0 2*10^2])
%     end
% end

%% GRs

GRs = gradient(real(log(blankeddata)), 10/60); % real

GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

strainGRs = zeros(nstrains, nmedia, nreplicates, ntimepoints);
conditioncounter = 2;
for i = 1:nconditions
    strainGRs(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = GRconditions{i+1}(:, :);
end
%%
% plot all replicates
figure(4)
hold on
for j = 1:nmedia
    for i = 1:nstrains
        subplot(nmedia, nstrains, (j-1)*nstrains+i)
        plot(timepoints, squeeze(strainGRs(i, j, :, :)), 'k');
        axis([0 24 -2 2])
    end
end

%% Averages

avdata = squeeze(mean(straindata, 3)); % strain plasmid media timepoint
avGRs = squeeze(mean(strainGRs, 3)); % 
avGFP = squeeze(mean(GFPstraindata, 3));
avBFP = squeeze(mean(BFPstraindata, 3));

stderrdata = squeeze(std(straindata, 0, 3)) / sqrt(nreplicates);
stderrGRs = squeeze(std(strainGRs, 0, 3)) / sqrt(nreplicates); % 0 is default weight
stderrGFP = squeeze(std(GFPstraindata, 0, 3)) / sqrt(nreplicates);
stderrBFP = squeeze(std(BFPstraindata, 0, 3)) / sqrt(nreplicates);
%%
% % plot all
% figure(5)
% hold on
% for j = 1:nmedia
%     for i = 1:nstrains
%         subplot(nmedia, nstrains, (j-1)*nstrains+i)
%         hold on
%         plot(timepoints, squeeze(avdata(i, j, :)), 'Color', 'k', 'LineWidth', 0.1);
%         errorbar(timepoints, squeeze(avdata(i, j, :)), squeeze(stderrdata(i, j, :)), 'Color', [0.5 0.5 0.5], 'CapSize', 0);
%         axis([0 24 0 1.5])
%     end
% end
% 
% figure(6)
% hold on
% for j = 1:nmedia
%     for i = 1:nstrains
%         subplot(nmedia, nstrains, (j-1)*nstrains+i)
%         hold on
%         plot(timepoints, squeeze(avGRs(i, j, :)), 'Color', 'k', 'LineWidth', 0.1);
%         errorbar(timepoints, squeeze(avGRs(i, j, :)), squeeze(stderrGRs(i, j, :)), 'Color', [0.5 0.5 0.5], 'CapSize', 0);
%         axis([0 24 -2 2])
%     end
% end


%% Figure 4A plots

cmap = cmocean('amp', 100);

straincolors = zeros(nstrains-1, 3);
straincolors(2:end, :) = cmap([75 10 90 25],:);

strainorder = [2 3 5 4 6];

% LB
figure(7)
hold on
handles = cell(nstrains-1, 1);
%for i = 2:nstrains
for i = strainorder
    plot(timepoints, squeeze(avdata(i, 1, :)), 'Color', straincolors(i-1, :), 'LineWidth', 7);
end
% for i = 2:nstrains
%    errorbar(timepoints, squeeze(avdata(i, 1, :)), squeeze(stderrdata(i, 1, :)), 'Color', straincolors(i-1, :) + 0.2*(ones(1, 3) - straincolors(i-1, :)), 'CapSize', 0);
% end

%xlabel("Hours");
%ylabel("OD");
xlim([0 24]);
ylim([0 1.5]);
%legend(strainnames(2:end), 'Location', 'eastoutside');
%legend(strainnames(strainorder), 'Location', 'eastoutside');
%title(medianames(1));
set(gca,'FontSize',30)
axis square
hold off

% LB with carbenicillin
figure(8)
hold on
%plot(timepoints, squeeze(avdata(2, 1, :)), 'Color', straincolors(1, :), 'LineWidth', 7, 'LineStyle', ':');
for i = strainorder
    plot(timepoints, squeeze(avdata(i, 2, :)), 'Color', straincolors(i-1, :), 'LineWidth', 7);
end
% errorbar(timepoints, squeeze(avdata(2, 1, :)), squeeze(stderrdata(1, 1, :)), 'Color', straincolors(1, :) + 0.2*(ones(1, 3) - straincolors(1, :)), 'CapSize', 0);
% for i = 2:nstrains
%     errorbar(timepoints, squeeze(avdata(i, 2, :)), squeeze(stderrdata(i, 2, :)), 'Color', straincolors(j, :) + 0.2*(ones(1, 3) - straincolors(j, :)), 'CapSize', 0);
% end
        
%xlabel("Hours");
%ylabel("OD");
xlim([0 24]);
ylim([0 1.5]);
%legend(strainnames(strainorder), 'Location', 'eastoutside');
%legend('boxoff');
%title(medianames(2));
set(gca,'FontSize',30)
axis square
hold off
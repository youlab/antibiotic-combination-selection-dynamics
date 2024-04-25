close all
clear

% set constants
nwells = 384;
height = 16;
width = 24;
ntimepoints = 145;
nreplicates = 2;
nmedia = 3;
nstrains = 49;

% read in raw data from excel
datafile = '2-20-2024 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D73:NW217')';
fluordata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D221:NW365')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B73:B217');

layoutfile = '2-20-2024 Plate Layout.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'E1:AB16');
layoutlong = reshape(platelayout', [nwells 1]);
conditioninfo = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'A2:C148');

medianames = ["LB" "CM 2 \mug/mL" "AMX 2 \mug/mL"];
strainnames = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'AE2:AE50', 'OutputType', 'string');

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
cleanblankdata = blankdata(max(blankdata, [], 2) < 0.1, :);
blanks = mean(cleanblankdata);
blankeddata = datasmoothed - blanks;
blankeddata(blankeddata<0) = 0;

% process fluor data
fdatasmooth = movmean(fluordata, 5, 2);
fblankdata = fdatasmooth(layoutlong==0, :);
cleanfblankdata = fblankdata(blankdata(:, 145) < 0.1, :);
fblanks = mean(cleanfblankdata);
fblankeddata = fdatasmooth - fblanks;
fblankeddata(fblankeddata<0) = 0;


%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

fconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    fconditionsdata{i+1} = fblankeddata(layoutlong==i, :);
end

straindata = nan(nstrains, nmedia, nreplicates, ntimepoints);
fstraindata = straindata;

for i = 1:nconditions
    straindata(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = conditionsdata{i+1}(:, :);
    fstraindata(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = fconditionsdata{i+1}(:, :);
end
%% GRs

GRs = gradient(real(log(blankeddata)), 10/60); % real

GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

strainGRs = nan(nstrains, nmedia, nreplicates, ntimepoints);
for i = 1:nconditions
    strainGRs(conditioninfo(i, 2), conditioninfo(i, 3), :, :) = GRconditions{i+1}(:, :);
end

%% averages

avdata = squeeze(mean(straindata, 3, 'omitnan')); % strain media timepoint
avGRs = squeeze(mean(strainGRs, 3, 'omitnan')); % 
avfluor = squeeze(mean(fstraindata, 3, 'omitnan'));

stderrdata = squeeze(std(straindata, 0, 3)) / sqrt(nreplicates); 
stderrGRs = squeeze(std(strainGRs, 0, 3)) / sqrt(nreplicates); % 0 is default weight
stderrfluor = squeeze(std(fstraindata, 0, 3)) / sqrt(nreplicates);

%% LB only

cmap = crameri('roma', nstrains);

figure(1)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avdata(i, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avdata(41, 1, :)), 'k', 'LineWidth', 5); %DA28102
xlabel("Hours");
ylabel("OD");
xlim([0 24]);
ylim([0 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%GRs
figure(2)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avGRs(i, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avGRs(41, 1, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("GR (1/h)");
xlim([0 10]);
ylim([-0.5 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%fluor
figure(3)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avfluor(i, 1, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avfluor(41, 1, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("mCherry");
xlim([0 24]);
ylim([0 70000]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off
%% LB + cm 2

figure(4)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avdata(i, 2, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avdata(31, 2, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("OD");
xlim([0 24]);
ylim([0 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%GRs
figure(5)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avGRs(i, 2, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avGRs(41, 2, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("GR (1/h)");
xlim([0 10]);
ylim([-0.5 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%fluor
figure(6)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avfluor(i, 2, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avfluor(41, 2, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("mCherry");
xlim([0 24]);
ylim([0 70000]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off
%% LB + amx 2

figure(7)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avdata(i, 3, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avdata(41, 3, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("OD");
xlim([0 24]);
ylim([0 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%GRs
figure(8)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avGRs(i, 3, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avGRs(41, 3, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("GR (1/h)");
xlim([0 10]);
ylim([-0.5 2]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%fluor
figure(9)
hold on
for i = 1:nstrains
    plot(timepoints, squeeze(avfluor(i, 3, :)), 'Color', cmap(i, :), 'LineWidth', 5)        
end
plot(timepoints, squeeze(avfluor(41, 3, :)), 'k', 'LineWidth', 5);
xlabel("Hours");
ylabel("mCherry");
xlim([0 24]);
ylim([0 70000]);
legend(strainnames, 'Location', 'eastoutside', 'NumColumns', 3);
set(gca,'FontSize',20)
axis square
hold off

%% 
cmap = crameri('lajolla', nmedia+4);

figure(10)
tiledlayout(7, 7, 'TileSpacing', 'tight', 'Padding', 'tight')
for i = 1:nstrains
    nexttile
    hold on
    for j = 1:nmedia
        plot(timepoints, squeeze(avdata(i, j, :)), 'Color', cmap(j+3, :), 'LineWidth', 3)
    end
    %xlabel("Hours");
    %ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    %title(strainnames(i))
    text(1, 1.8, strainnames(i))
    %legend(medianames, 'Location', 'eastoutside');
    %set(gca,'FontSize',20)
    axis square
    hold off
    %set(fig, 'WindowState', 'maximized')
    %saveas(fig, strainnames(i), 'jpeg')
end

% replot for 30 strains only
% wrong names (correct wells)
%strains30 = ["K1", "K2", "K9", "K12", "K19", "K25", "K29", "K30", "K31", "K32", "K33", "K45", "K47", "K49", "K51", "K53", "K54", "K58", "K64", "K65", "K67", "K68", "K70", "K78", "K81", "K82", "K84", "K85", "K89", "K96"];
strains30 = ["K1", "K2", "K4", "K5", "K7", "K9", "K11", "K12", "K13", "K20", "K23", "K27", "K31", "K33", "K38", "K45", "K47", "K52", "K54", "K55", "K59", "K64", "K67", "K69", "K70", "K76", "K78", "K88", "K92", "K96"];
strain30indices = matches(strainnames, strains30);
strainindices = 1:nstrains;

figure(11)
tiledlayout(5, 6, 'TileSpacing', 'tight', 'Padding', 'tight')
for i = strainindices(strain30indices)
    nexttile
    hold on
    for j = 1:nmedia
        plot(timepoints, squeeze(avdata(i, j, :)), 'Color', cmap(j+3, :), 'LineWidth', 3)
    end
    %xlabel("Hours");
    %ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    %title(strainnames(i))
    text(1, 1.8, strainnames(i))
    %legend(medianames, 'Location', 'eastoutside');
    %set(gca,'FontSize',20)
    axis square
    hold off
    %set(fig, 'WindowState', 'maximized')
    %saveas(fig, strainnames(i), 'jpeg')
end
%%
%legend / spot testing
figure(12)
hold on
for j = 1:nmedia
    plot(timepoints, squeeze(avdata(5, j, :)), 'Color', cmap(j+3, :), 'LineWidth', 3)
end
%xlabel("Hours");
%ylabel("OD");
xlim([0 24]);
ylim([0 2]);
%title(strainnames(i))
text(1, 1.8, strainnames(i))
legend(medianames, 'Location', 'eastoutside');
legend boxoff;
set(gca,'FontSize',20)
axis square
hold off
%set(fig, 'WindowState', 'maximized')
%saveas(fig, strainnames(i), 'jpeg')

%% max GRs

maxCMGRs = max(squeeze(avGRs(:, 2, :)), [], 2);


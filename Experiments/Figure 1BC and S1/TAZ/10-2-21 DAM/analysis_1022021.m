close all
clear

% set constants
nwells = 384;
height = 16;
width = 24;
ntimepoints = 145;
antibiotics = [0 0.5 1 2 4 8 16 32 64 128];
inhibitors = antibiotics;
nconditions = length(antibiotics)*length(inhibitors);
nreplicates = floor(nwells / nconditions);

strainname = 'DA28102 Low Copy BlaM';

% read in raw data from excel
datafile = '10-2-2021 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D92:NW236')';
gfpdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D240:NW384')';
bfpdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D388:NW532')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B92:B236');

layoutfile = '10-2-2021 Protocol.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'E1:AB16');
layoutlong = reshape(platelayout', [nwells 1]);

timepoints = timedata ./ 3600;

%% process data

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% smooth
datasmoothed = movmean(rawdata, 5, 2);

% blanking smoothed data
blankdata = datasmoothed(layoutlong == 0, :);
% threshold the blanks to remove contamination
cleanblankdata = blankdata(blankdata(:, 145) < 0.1, :);
blanks = mean(cleanblankdata);
blankeddata = datasmoothed - blanks;

% process gfp data
gfpdatasmooth = movmean(gfpdata, 5, 2);
gfpblankdata = gfpdatasmooth(layoutlong==0, :);
cleangfpblankdata = gfpblankdata(blankdata(:, 145) < 0.1, :);
gfpblanks = mean(cleangfpblankdata);
gfpblankeddata = gfpdatasmooth - gfpblanks;

% process bfp data
bfpdatasmooth = movmean(bfpdata, 5, 2);
bfpblankdata = bfpdatasmooth(layoutlong==0, :);
cleanbfpblankdata = bfpblankdata(blankdata(:, 145) < 0.1, :);
bfpblanks = mean(cleanbfpblankdata);
bfpblankeddata = bfpdatasmooth - bfpblanks;

%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

gfpconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    gfpconditionsdata{i+1} = gfpblankeddata(layoutlong==i, :);
end

bfpconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    bfpconditionsdata{i+1} = bfpblankeddata(layoutlong==i, :);
end

concentrationdata = zeros(length(antibiotics), length(inhibitors), nreplicates, ntimepoints);
gfpconcentrationdata = concentrationdata;
bfpconcentrationdata = concentrationdata;

conditioncounter = 2; % start after blanks
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
            concentrationdata(i, j, :, :) = conditionsdata{conditioncounter}(:, :);
            gfpconcentrationdata(i, j, :, :) = gfpconditionsdata{conditioncounter}(:, :);
            bfpconcentrationdata(i, j, :, :) = bfpconditionsdata{conditioncounter}(:, :);
            conditioncounter = conditioncounter + 1;
    end
end

%%
% plot all replicates
figure(1)
hold on
counter = 1;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
            subplot(10, 10, counter)
            plot(timepoints, squeeze(concentrationdata(i, j, :, :)), 'k');
            counter = counter+1;
            axis([0 24 0 1.5])
    end
end

figure(2)
hold on
counter = 1;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
        subplot(10, 10, counter)
        plot(timepoints, squeeze(gfpconcentrationdata(i, j, :, :)), 'k');
        counter = counter+1;
        axis([0 24 0 4*10^3])
    end
end

figure(3)
hold on
counter = 1;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
        subplot(10, 10, counter)
        plot(timepoints, squeeze(bfpconcentrationdata(i, j, :, :)), 'k');
        counter = counter+1;
        axis([0 24 0 50])
    end
end

%% GRs

GRs = gradient(log(blankeddata), 10/60); % real

GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

concentrationGRs = zeros(length(antibiotics), length(inhibitors), nreplicates, ntimepoints);
conditioncounter = 2;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
        concentrationGRs(i, j, :, :) = GRconditions{conditioncounter}(:, :);
        conditioncounter = conditioncounter + 1;
    end
end
%%
% plot all replicates
figure(4)
hold on
counter = 1;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
        subplot(10, 10, counter)
        plot(timepoints, squeeze(concentrationGRs(i, j, :, :)), 'k');
        counter = counter+1;
        axis([0 24 -2 2])
    end
end

%% Averages

avdata = squeeze(mean(concentrationdata, 3)); % [a] [i] timepoint
avGRs = squeeze(mean(concentrationGRs, 3)); % 
avgfp = squeeze(mean(gfpconcentrationdata, 3));
avbfp = squeeze(mean(bfpconcentrationdata, 3));

stderrdata = squeeze(std(concentrationdata, 0, 3)) / sqrt(nreplicates);
stderrGRs = squeeze(std(concentrationGRs, 0, 3)) / sqrt(nreplicates); % 0 is default weight
stderrgfp = squeeze(std(gfpconcentrationdata, 0, 3)) / sqrt(nreplicates);
stderrbfp = squeeze(std(bfpconcentrationdata, 0, 3)) / sqrt(nreplicates);

%% plot condition timecourses with average and stderr

% plot all replicates
figure(5)
hold on
counter = 1;
for i = 1:length(antibiotics)
    for j = 1:length(inhibitors)
            subplot(10, 10, counter)
            hold on
            plot(timepoints, squeeze(avdata(i, j, :)), 'k', 'LineWidth', 0.2);
            errorbar(timepoints, squeeze(avdata(i, j, :)), squeeze(stderrdata(i, j, :)), 'Color', [0.8 0.8 0.8], 'CapSize', 0);
            counter = counter+1;
            axis([0 24 0 1.5])
            hold off
    end
end


%% heatmap
for hour = [4 8 16 24]
    
    finalOD = squeeze(avdata(:, :, find(abs(timepoints - hour) < 0.05, 1, 'first')));
    finalgfp = squeeze(avgfp(:, :, find(abs(timepoints - hour) < 0.05, 1, 'first')));

    figure(100+hour)
    imagesc(finalOD)
    %colorbar
    caxis([0 2])
    colormap(cmocean('tempo'))
    ax = gca;
    ax.YDir = 'normal';
    ax.XTick = linspace(1, 10, length(inhibitors)/3+1);
    ax.YTick = linspace(1, 10, length(antibiotics)/3+1);
    ax.XTickLabel = string(inhibitors(1:3:end));
    ax.YTickLabel = string(antibiotics(1:3:end)); 
    %xlabel("Tazobactam (\mug/mL)")
    %ylabel("Amoxicillin (\mug/mL)")
    %title(strainname)% + " " + hour +"h Cell Density")
    set(gca, 'fontsize', 30)
    axis square
    hold off
    set(gcf, 'position', [500 300 700 500])
end


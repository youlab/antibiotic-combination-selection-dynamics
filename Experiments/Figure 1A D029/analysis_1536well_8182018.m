%% analyze 99-repeat 1536-well plate on Spark
% DICON 029

% set constants
nwells = 1536;
height = 32;
width = 48;
ntimepoints = 121;
nconditions = 256; 
nreplicates = 5;
nantibiotic = 16;
ninhibitor = 16;
nblanks = 128;
ncontrols = 128;
datastartsrow = 5;
growththreshold = 0.3; %%%%%%% change this as needed
rawinterval = height + datastartsrow;

% read in raw data from excels
datafile = '8-18-2018 D029 Test Data';
rawdata = xlsread(datafile, 1, 'B49:AW4524');
protocolfile = '8-18-2018 D029 Test Protocol.xlsx';

conditions_concentrations = xlsread(protocolfile, 2, 'A:C');
platelayout = xlsread(protocolfile, 2, 'F1:BA32');

%% Extract Data

% extract timepoints
timedata = rawdata(1:rawinterval:end, 1);
timepoints = timedata ./ 3600;

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% save data by well
wellsdata = zeros(nwells, ntimepoints);
for i = 1:nwells
    for j = 1:ntimepoints
        [row, col] = find(wellindices == i);
        wellsdata(i, j) = rawdata(datastartsrow + (row - 1) + (rawinterval * (j - 1)), col);
    end
end

% save data by condition
conditionsdata = zeros(nconditions, nreplicates, ntimepoints);
for i = 1:nconditions
    [conditionrow, conditioncol] = find(platelayout == i);
    conditionindices = zeros(nreplicates, 1);
    for j = 1:nreplicates
        conditionindices(j) = wellindices(conditionrow(j), conditioncol(j));
    end
    conditionsdata(i, :, :) = wellsdata(conditionindices, :);
end

% collect blanks
blanksdata = wellsdata(wellindices(platelayout == -1), :);
blanksfinal = squeeze(blanksdata(:, ntimepoints)); % final OD of blanks
nblankscontam = sum(blanksfinal > growththreshold); % contaminated
blanksaverage = mean(blanksdata(~(blanksfinal > growththreshold), :)); % average of uncontaminated wells
blanksaveragesmoothed = movmean(blanksaverage, 5, 2);
% plot(timepoints, blanksaveragesmoothed); 

% collect controls
controlsdata = wellsdata(wellindices(platelayout == 0), :);
controlsfinal = squeeze(controlsdata(:, ntimepoints));
controlsaverage = mean(controlsdata(~(controlsfinal < growththreshold), :)); % average of growing wells
controlsaveragesmoothed = movmean(controlsaverage, 5, 2);
% plot(timepoints, controlsaveragesmoothed);

%% Process Data

wellsdatasmoothed = movmean(wellsdata, 5, 2);
wellsdatafinal = wellsdatasmoothed - repmat(blanksaveragesmoothed, nwells, 1);

% save smoothed and blanked by condition
conditionsdatafinal = zeros(nconditions, nreplicates, ntimepoints);
for i = 1:nconditions
    [conditionrow, conditioncol] = find(platelayout == i);
    conditionindices = zeros(nreplicates, 1);
    for j = 1:nreplicates
        conditionindices(j) = wellindices(conditionrow(j), conditioncol(j));
    end
    conditionsdatafinal(i, :, :) = wellsdatafinal(conditionindices, :);
end

% average conditions
conditionsaverage = squeeze(mean(conditionsdata, 2));
conditionsaveragefinal = squeeze(mean(conditionsdatafinal, 2));

%% extract heatmap from conditions

finalODconditions = conditionsaveragefinal(:, ntimepoints);
antibiotic = conditions_concentrations((1 + ninhibitor*((1:nantibiotic) - 1)), 2);
inhibitor = conditions_concentrations(1:ninhibitor, 3);
finalODmap = reshape(finalODconditions, [ninhibitor nantibiotic])';
%%
figure(1)
imagesc([min(inhibitor) max(inhibitor)], [max(antibiotic) min(antibiotic)], flipud(finalODmap));
%colorbar
caxis([0 2]) 
colormap(cmocean('tempo'))
ax = gca;
axis square
ax.YDir = 'normal';
%labels
%title('DICON 029')
%xlabel('clavulanic acid (\mug/mL)')    
%ylabel('cefotaxime (\mug/mL)')
ax.XTick = linspace(min(inhibitor), max(inhibitor), ninhibitor/5 + 1);
ax.YTick = linspace(min(antibiotic), max(antibiotic), nantibiotic/5 + 1);
indices = 1:(ninhibitor-1)/3:16;
ax.XTickLabel = compose('%-.2g', inhibitor(indices));
ax.YTickLabel = string(antibiotic(indices));
set(gca, 'fontsize', 30')
set(gcf, 'position', [500 300 700 500])
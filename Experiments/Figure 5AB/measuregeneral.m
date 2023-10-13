% measure.m calculates many different parameters from an arbitrary number
% of growth curves.

function [growthRates, crashData, peakData, RTData] = measuregeneral(growthdatafull, iterationTime)

% get time and plate information
ntimepoints = size(growthdatafull, 2);
nwells = size(growthdatafull, 1);
timepoints = (0:ntimepoints-1)*iterationTime;

% remove noisy values
growthdatafull(growthdatafull < .005) = NaN;

% log derivatives
growthRates = gradient(log(growthdatafull), 10/60);

% Find maximum recovery rate, the corresponding OD, and the time to that point 
peakData = zeros(nwells, 3);
skipFirstPoints = 24; % this is to avoid detecting a peak in the noisy area (where values are still small)
[peakData(:, 1), iGRPeak] = max(growthRates(:,skipFirstPoints+1:end), [], 2);
iGRPeak = iGRPeak + skipFirstPoints;
peakData(:, 2) = timepoints(iGRPeak);
idx = sub2ind(size(growthdatafull), 1:nwells, iGRPeak');
peakData(:, 3) = growthdatafull(idx); 

% Find maximum crash rate, the corresponding OD, and the time to that point
crashData = zeros(nwells, 3);
skipFirstPoints = 10;
[crashData(:, 1), iGRCrash] = min(growthRates(:,skipFirstPoints+1:43), [], 2);
iGRCrash=iGRCrash+skipFirstPoints;
crashData(:, 2) = timepoints(iGRCrash);
idx = sub2ind(size(growthdatafull), 1:nwells, iGRCrash');
crashData(:, 3) = growthdatafull(idx); 

% Find recovery time, the corresponding OD, and the time to that point
RTData = zeros(nwells, 2);
[maxOD, ~]=max(growthdatafull, [], 2);
threshold=0.5.*maxOD;
[~, iRT] = max(growthdatafull >= threshold, [], 2); % finds index of first point over the threshold
RTData(:, 1) = timepoints(iRT); 
idx = sub2ind(size(growthdatafull), 1:nwells, iRT');
RTData(:, 2) = growthdatafull(idx);
RTData((RTData(:, 1) == 0), 1) = Inf; % if it never recovered
RTData((RTData(:, 1) == 0), 2) = NaN;
end 
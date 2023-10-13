close all
clear

% set constants
nwells = 384;
height = 16;
width = 24;
ntimepoints = 145;
nconditions = 72;

% read in raw data from excel
datafile = '1-13-2020 2004 2055 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D54:NW198')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B54:B198');

layoutfile = '1-13-2020 2004 2055 Protocol.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'G1:AD16');
layoutlong = reshape(platelayout', [384 1]);

timepoints = timedata ./ 3600;

%% process data

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% smooth
datasmoothed = movmean(rawdata, 5, 2);

% get blanks
blanks = mean(datasmoothed(layoutlong == 0, :));
blankeddata = datasmoothed - blanks;

% blanks = mean(rawdata(layoutlong == 0, :));
% blankeddata = rawdata - blanks;

%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

% testcondition1 = conditionsdata{71}; %here, the thing in parentheses is the condition - 1
% testcondition2 = conditionsdata{1};

% % plot a specific condition
% figure(6);
% hold on
% for i = 1:size(testcondition1, 1)
%     plot(timepoints, testcondition1(i, :), 'r');
% end
% for i = 1:size(testcondition2, 1)
%     plot(timepoints, testcondition2(i, :), 'g');
% end
% xlabel("Hours");
% ylabel("OD");
% ylim([0 2]);
% xlim([0 24]);
% set(gca,'FontSize',20)
% hold off;

% % plot all conditions
% figure(21)
% hold on
% for i = 1:nconditions+1
%     linecolor = rand(1, 3);
%     for j = 1:size(conditionsdata{i}, 1)
%         plot(timepoints, conditionsdata{i}(j, :), 'Color', linecolor)
%     end
% end
% hold off

% noantibiotic = zeros((nconditions / 3)*4, ntimepoints);
% for i =  1:(nconditions / 3)
%     noantibiotic(4*i - 3:4*i, :) = conditionsdata{3*i - 1};
% end
% figure(2322232)
% plot(timepoints, noantibiotic)
%% plot conditions

% make isolate colors
colors = [215 48 39;
        244 109 67;
        253 174 97;
        254 224 144;
        224 243 248;
        171 217 233;
        116 173 209;
        69 117 180];
% colors = [140 81 10;
%           191 129 45;
%           223 194 125;
%           246 232 195;
%           199 234 229;
%           128 205 193;
%           53 151 143;
%           1 102 94];
colors = colors ./ 255;

% plot all no-antibiotic conditions (carolyn's)
figure(2836592)
hold on
for i = 2:nconditions+1
    if mod(i, 3) == 2        
        linecolor = colors(floor((i - 2) / 9) + 1, :); 
        for j = 1:size(conditionsdata{i}, 1)
            plot(timepoints, conditionsdata{i}(j, :), 'Color', linecolor)
        end
    end
end
xlabel("Hours");
ylabel("OD");
ylim([0 2]);
xlim([0 24]);
set(gca,'FontSize',20)
title("No Antibiotic")
hold off

% plot all cefotaxime only conditions
figure(283593222)
hold on
for i = 2:nconditions+1
    if mod(i, 3) == 0    
        linecolor = colors(floor((i - 2) / 9) + 1, :); 
        for j = 1:size(conditionsdata{i}, 1)
            plot(timepoints, conditionsdata{i}(j, :), 'Color', linecolor)
        end
    end
end
xlabel("Hours");
ylabel("OD");
ylim([0 2]);
xlim([0 24]);
set(gca,'FontSize',20)
title("50 \mug/mL CTX")
hold off

% plot all cefo + clav conditions
figure(222325225)
hold on
for i = 2:nconditions+1
    if mod(i, 3) == 1    
        linecolor = colors(floor((i - 2) / 9) + 1, :); 
        for j = 1:size(conditionsdata{i}, 1)
            plot(timepoints, conditionsdata{i}(j, :), 'Color', linecolor)
        end
    end
end
xlabel("Hours");
ylabel("OD");
ylim([0 2]);
xlim([0 24]);
set(gca,'FontSize',20)
title("50 \mug/mL CTX + 0.5 \mug/mL CLA")
hold off
%% final heatmap
finalOD = blankeddata(:, 145);
finalplate = reshape(finalOD, [width height])';
figure(10)
hold on
imagesc(flipud(finalplate))
colorbar
axis tight;
hold off
%% heatmap with contamination

count = 0;
contamfinalOD = finalOD;
for i = 1:nwells
    if finalOD(i) > 0.08
        count = count + 1;
    end
    if (finalOD(i) > 0.08 && layoutlong(i) == 0)
        contamfinalOD(i) = -1;
        disp(i);
    end
end

finalplatecontam = reshape(contamfinalOD, [width height])';
figure(11)
hold on
imagesc(flipud(finalplatecontam))
colorbar
axis tight;
hold off

%% GRs

GRs = gradient(blankeddata, 10/60); %%% CHANGE BASED ON RAW VS SMOOTHED
GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

testcondition1 = GRconditions{1};
testcondition2 = GRconditions{71};

% plot a specific condition
figure(38692);
hold on
for i = 1:size(testcondition1, 1)
    plot(timepoints, testcondition1(i, :), 'r');
end
for i = 1:size(testcondition2, 1)
    plot(timepoints, testcondition2(i, :), 'g');
end
xlabel("Hours");
ylabel("d(OD)/dt");
xlim([0 24]);
set(gca,'FontSize',20)
hold off;
%%
% figure(295638)
% hold on
% for i = 1:size(GRconditions{65}, 1)
%     plot(timepoints, GRconditions{65}(i, :), 'r');
% end
% for i = 1:size(GRconditions{68}, 1)
%     plot(timepoints, GRconditions{68}(i, :), 'g');
% end
% for i = 1:size(GRconditions{71}, 1)
%     plot(timepoints, GRconditions{71}(i, :), 'b');
% end
% xlabel("Hours");
% ylabel("d(OD)/dt");
% xlim([0 24]);
% ylim([0 0.4]);
% set(gca,'FontSize',20)
% title('2050')
% hold off;

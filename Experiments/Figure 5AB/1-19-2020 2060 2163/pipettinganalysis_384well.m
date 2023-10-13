close all
clear

% set constants
nwells = 384;
height = 16;
width = 24;
ntimepoints = 145;
nconditions = 72;
nisolates = 8;
ntreatments = 3;
ntechreps = 4;
nbioreps = 3;

isolatenames = ["2060", "2091", "2094", "2099", "2137", "2140", "2148", "2163"];

% read in raw data from excel
datafile = '1-19-2020 2060 2163 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D54:NW198')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B54:B198');

layoutfile = '1-19-2020 2060 2163 Protocol.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'G1:AD16');
layoutlong = reshape(platelayout', [384 1]);

timepoints = timedata ./ 3600;

%% process data

% numerically index the plate layout
wellindices = reshape((1:nwells), [width height])';

% smooth
datasmoothed = movmean(rawdata, 5, 2);

% blanking smoothed data
blanks = mean(datasmoothed(layoutlong == 0, :));
blankeddata = datasmoothed - blanks;

% blanking raw data
rblanks = mean(rawdata(layoutlong == 0, :));
rblankeddata = rawdata - rblanks;

%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

rconditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    rconditionsdata{i+1} = rblankeddata(layoutlong == i, :);
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

% carolyn's version
noantibiotic = zeros((nconditions / 3)*4, ntimepoints);
for i =  1:(nconditions / 3)
    noantibiotic(4*i - 3:4*i, :) = rconditionsdata{3*i - 1};
end
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

% plot all no-antibiotic conditions
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

% plot all amoxicillin only conditions
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
title("50 \mug/mL AMX")
hold off

% plot all amox + clav conditions
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
title("50 \mug/mL AMX + 25 \mug/mL CLA")
hold off

%% Plot individual isolates

isolatedata = zeros(nisolates, ntreatments, nbioreps*ntechreps, ntimepoints);

for i = 1:nisolates
    noantibiotic = cat(1, conditionsdata{9*i - 7}, conditionsdata{9*i - 4}, conditionsdata{9*i - 1});
    amx = cat(1, conditionsdata{9*i - 6}, conditionsdata{9*i - 3}, conditionsdata{9*i});
    amxcla = cat(1, conditionsdata{9*i - 5}, conditionsdata{9*i - 2}, conditionsdata{9*i + 1});
    
    isolatedata(i, 1, :, :) = noantibiotic;
    isolatedata(i, 2, :, :) = amx;
    isolatedata(i, 3, :, :) = amxcla;
    
    figure(str2double(isolatenames(i))+10000)
    hold on
    
    for j = 1:size(noantibiotic, 1)
        set1 = plot(timepoints, noantibiotic(j, :), 'Color', colors(4, :), 'LineWidth', 2);
    end
    for j = 1:size(amx, 1)
        set2 = plot(timepoints, amx(j, :), 'Color', colors(1, :), 'LineWidth', 2);
    end
    for j = 1:size(amxcla, 1)
        set3 = plot(timepoints, amxcla(j, :), 'Color', colors(8, :), 'LineWidth', 2);
    end
    
    xlabel("Hours");
    ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    title(isolatenames(i))
    legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off;
end
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

dGRs = gradient(blankeddata, 10/60); % Carolyn's
GRs = gradient(log(blankeddata), 10/60); % real

dGRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    dGRconditions{i+1} = dGRs(layoutlong == i, :);
end
GRconditions = cell(nconditions+1, 1);
for i = 0:nconditions
    GRconditions{i+1} = GRs(layoutlong == i, :);
end

% testcondition1 = dGRconditions{1};
% testcondition2 = dGRconditions{71};

% % plot a specific condition
% figure(38692);
% hold on
% for i = 1:size(testcondition1, 1)
%     plot(timepoints, testcondition1(i, :), 'r');
% end
% for i = 1:size(testcondition2, 1)
%     plot(timepoints, testcondition2(i, :), 'g');
% end
% xlabel("Hours");
% ylabel("d(OD)/dt");
% xlim([0 24]);
% set(gca,'FontSize',20)
% hold off;
%% Carolyn plotting

% for i = 1:nisolates
%     bioreplicate1 = dGRconditions{9*i - 7};
%     bioreplicate2 = dGRconditions{9*i - 4};
%     bioreplicate3 = dGRconditions{9*i - 1};
%     
%     figure(str2double(isolatenames(i)))
%     hold on
%     
%     for j = 1:size(bioreplicate1, 1)
%         plot(timepoints, bioreplicate1(j, :), 'r');
%     end
%     for j = 1:size(bioreplicate2, 1)
%         plot(timepoints, bioreplicate2(j, :), 'g');
%     end
%     for j = 1:size(bioreplicate3, 1)
%         plot(timepoints, bioreplicate3(j, :), 'b');
%     end
%     
%     xlabel("Hours");
%     ylabel("d(OD)/dt");
%     xlim([0 24]);
%     ylim([-0.3 0.7]);
%     set(gca,'FontSize',20)
%     title(isolatenames(i))
%     hold off;
% end

%% my GR plotting

isolateGRs = zeros(nisolates, ntreatments, nbioreps*ntechreps, ntimepoints);

for i = 1:nisolates
    noantibiotic = cat(1, GRconditions{9*i - 7}, GRconditions{9*i - 4}, GRconditions{9*i - 1});
    amx = cat(1, GRconditions{9*i - 6}, GRconditions{9*i - 3}, GRconditions{9*i});
    amxcla = cat(1, GRconditions{9*i - 5}, GRconditions{9*i - 2}, GRconditions{9*i + 1});
    
    isolateGRs(i, 1, :, :) = noantibiotic;
    isolateGRs(i, 2, :, :) = amx;
    isolateGRs(i, 3, :, :) = amxcla;
    
    figure(str2double(isolatenames(i))+20000)
    hold on
    
    for j = 1:size(noantibiotic, 1)
        set1 = plot(timepoints, noantibiotic(j, :), 'Color', colors(4, :), 'LineWidth', 2);
    end
    for j = 1:size(amx, 1)
        set2 = plot(timepoints, amx(j, :), 'Color', colors(1, :), 'LineWidth', 2);
    end
    for j = 1:size(amxcla, 1)
        set3 = plot(timepoints, amxcla(j, :), 'Color', colors(8, :), 'LineWidth', 2);
    end
    
    xlabel("Hours");
    ylabel("GR");
    xlim([0 24]);
    ylim([-0.5 1.2]);
    title(isolatenames(i))
    legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off;
end
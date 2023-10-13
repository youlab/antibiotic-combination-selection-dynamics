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
ntotalreps = ntechreps * nbioreps;

isolatenames = ["DICON185", "DICON193", "DICON195", "DICON205", "DICON208", "DICON210", "DICON211", "blank"];

% read in raw data from excel
datafile = '10-12-2020 DICON185 DICON211 Data.xlsx';
rawdata = readmatrix(datafile, 'Sheet', 1, 'Range', 'D54:NW198')';
timedata = readmatrix(datafile, 'Sheet', 1, 'Range', 'B54:B198');

layoutfile = '10-12-2020 DICON185 DICON211 Protocol.xlsx';
platelayout = readmatrix(layoutfile, 'Sheet', 2, 'Range', 'G1:AD16');
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

% blanking raw data
rblankdata = rawdata(layoutlong == 0, :);
cleanrblankdata = rblankdata(rblankdata(:, 145) < 0.1, :);
rblanks = mean(cleanrblankdata);
rblankeddata = rawdata - rblanks;

%% pull out data for each condition

conditionsdata = cell(nconditions+1, 1);
for i = 0:nconditions
    conditionsdata{i+1} = blankeddata(layoutlong == i, :);
end

rconditionsdata = cell(nconditions+1, 1); % raw conditions data
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

% carolyn's version, no smoothing
carolyndata = zeros((nconditions / 3)*4, ntimepoints);
for i =  1:(nconditions / 3)
    carolyndata(4*i - 3:4*i, :) = rconditionsdata{3*i - 1};
end
cdatanames = repmat(isolatenames, ntotalreps, 1);
cdatanames = reshape(cdatanames, ntotalreps*nisolates, 1);

% yasa data
yasadata = zeros(nconditions*4, ntimepoints);
for i = 1:nconditions % conditions go: bio rep, treatment, tech rep
    yasadata(4*i-3:4*i, :) = rconditionsdata{i + 1};
end
ydatanames = repelem(isolatenames, ntotalreps*ntreatments)';
yconditioncodes = zeros(nconditions*4, 3);
yconditioncodes(:, 1) = repmat(repelem(1:nbioreps, ntechreps*ntreatments), 1, nisolates); % bio rep
yconditioncodes(:, 2) = repmat(repelem(1:ntreatments, ntechreps), 1, nisolates*nbioreps); % treatment
yconditioncodes(:, 3) = repmat(1:ntechreps, 1, nisolates*ntreatments*nbioreps); % tech rep

filename = strcat("blankeddata ", isolatenames(1), " ", isolatenames(length(isolatenames)), ".xlsx");
writematrix(ydatanames, filename, 'Sheet', 1, 'Range', 'A1:A288');
writematrix(yconditioncodes, filename, 'Sheet', 1, 'Range', 'B1:D288');
writematrix(yasadata, filename, 'Sheet', 1, 'Range', 'E1:ES288');
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
fig = figure(2836592);
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

saveas(fig, "no abx", "jpeg")

% plot all amoxicillin only conditions
fig = figure(283593222);
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

saveas(fig, "amx", "jpeg")

% plot all amox + clav conditions
fig = figure(222325225);
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

saveas(fig, "amx-cla", "jpeg")


%% Plot individual isolates

isolatedata = zeros(nisolates, ntreatments, ntotalreps, ntimepoints);

for i = 1:nisolates
    noantibiotic = cat(1, conditionsdata{9*i - 7}, conditionsdata{9*i - 4}, conditionsdata{9*i - 1});
    amx = cat(1, conditionsdata{9*i - 6}, conditionsdata{9*i - 3}, conditionsdata{9*i});
    amxcla = cat(1, conditionsdata{9*i - 5}, conditionsdata{9*i - 2}, conditionsdata{9*i + 1});
    
    isolatedata(i, 1, :, :) = noantibiotic;
    isolatedata(i, 2, :, :) = amx;
    isolatedata(i, 3, :, :) = amxcla;
    
    fig = figure(i+10000);
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
    
    saveas(fig, strcat(isolatenames(i), "OD"), "jpeg");
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

dGRs = gradient(rblankeddata, 10/60); % Carolyn's
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

for i = 1:nisolates
    bioreplicate1 = dGRconditions{9*i - 7};
    bioreplicate2 = dGRconditions{9*i - 4};
    bioreplicate3 = dGRconditions{9*i - 1};
    
    fig = figure(i);
    hold on
    
    for j = 1:size(bioreplicate1, 1)
        plot(timepoints, bioreplicate1(j, :), 'r');
    end
    for j = 1:size(bioreplicate2, 1)
        plot(timepoints, bioreplicate2(j, :), 'g');
    end
    for j = 1:size(bioreplicate3, 1)
        plot(timepoints, bioreplicate3(j, :), 'b');
    end
    
    xlabel("Hours");
    ylabel("d(OD)/dt");
    xlim([0 24]);
    ylim([-0.3 0.7]);
    set(gca,'FontSize',20)
    title(isolatenames(i))
    hold off;
    
    saveas(fig, strcat(isolatenames(i), "Derivative"), "jpeg");
end

%% my GR plotting

isolateGRs = zeros(nisolates, ntreatments, ntotalreps, ntimepoints);

for i = 1:nisolates
    noantibiotic = cat(1, GRconditions{9*i - 7}, GRconditions{9*i - 4}, GRconditions{9*i - 1});
    amx = cat(1, GRconditions{9*i - 6}, GRconditions{9*i - 3}, GRconditions{9*i});
    amxcla = cat(1, GRconditions{9*i - 5}, GRconditions{9*i - 2}, GRconditions{9*i + 1});
    
    isolateGRs(i, 1, :, :) = noantibiotic;
    isolateGRs(i, 2, :, :) = amx;
    isolateGRs(i, 3, :, :) = amxcla;
    
    fig = figure(i+20000);
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
    ylim([-0.5 1.5]);
    title(isolatenames(i))
    legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off;
    
    saveas(fig, strcat(isolatenames(i), "GR"), "jpeg");
end

%% Averages

avisoGRs = squeeze(mean(isolateGRs, 3)); % isolate x treatment (no, amx, amxcla) x timepoint
avisodata = squeeze(mean(isolatedata, 3)); % isolate x treatment (no, amx, amxcla) x timepoint

stdevisoGRs = squeeze(std(isolateGRs, 0, 3)) / sqrt(ntotalreps); % 0 is default weight
stdevisodata = squeeze(std(isolatedata, 0, 3)) / sqrt(ntotalreps);

%GRs
for i = 1:nisolates   
    fig = figure(i+30000);
    hold on
    
    set1 = plot(timepoints, squeeze(avisoGRs(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
    set2 = plot(timepoints, squeeze(avisoGRs(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
    set3 = plot(timepoints, squeeze(avisoGRs(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
    
    err1 = errorbar(timepoints, squeeze(avisoGRs(i, 1, :)), squeeze(stdevisoGRs(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
    err2 = errorbar(timepoints, squeeze(avisoGRs(i, 2, :)), squeeze(stdevisoGRs(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
    err3 = errorbar(timepoints, squeeze(avisoGRs(i, 3, :)), squeeze(stdevisoGRs(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
    
    xlabel("Hours");
    ylabel("GR");
    xlim([0 24]);
    ylim([-0.5 1.5]);
    title(isolatenames(i))
    legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off;
    
    saveas(fig, strcat(isolatenames(i), "GRAverage"), "jpeg");
end

% ODs
for i = 1:nisolates 
    fig = figure(i+40000);
    hold on
    
    set1 = plot(timepoints, squeeze(avisodata(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
    set2 = plot(timepoints, squeeze(avisodata(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
    set3 = plot(timepoints, squeeze(avisodata(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
    
    errorbar(timepoints, squeeze(avisodata(i, 1, :)), squeeze(stdevisodata(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
    errorbar(timepoints, squeeze(avisodata(i, 2, :)), squeeze(stdevisodata(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
    errorbar(timepoints, squeeze(avisodata(i, 3, :)), squeeze(stdevisodata(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
    
    xlabel("Hours");
    ylabel("OD");
    xlim([0 24]);
    ylim([0 2]);
    title(isolatenames(i))
    legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
    set(gca,'FontSize',20)
    axis square
    hold off;
    
    saveas(fig, strcat(isolatenames(i), "ODAverage"), "jpeg");
end


%% Measurements?

% % private benefit by max no antibiotic peak
 privatebenefit = zeros(nisolates, 5); % GR, OD, timepoint, corresponding abx GR, calculation

[privatebenefit(:, 1), privatebenefit(:, 3)] = max(squeeze(avisoGRs(:, 1, :)), [], 2);
for i = 1:nisolates
    privatebenefit(i, 2) = avisodata(i, 1, privatebenefit(i, 3));
    privatebenefit(i, 4) = avisoGRs(i, 2, privatebenefit(i, 3));
    privatebenefit(i, 5) = privatebenefit(i, 4) / privatebenefit(i, 1);
end
% 
% b = bar(privatebenefit(:, 5));
% b.FaceColor = 'flat';
% b.CData = colors;
% set(gca, 'XTickLabel', isolatenames);
% ylabel("Betamin by no-abx peak")
% set(gca, 'FontSize', 20);
% 
% % plot peak
% %GRs
% for i = 1:nisolates   
%     figure(str2double(isolatenames(i))+50000)
%     hold on
%     
%     set1 = plot(timepoints, squeeze(avisoGRs(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
%     set2 = plot(timepoints, squeeze(avisoGRs(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
%     set3 = plot(timepoints, squeeze(avisoGRs(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
%     
%     err1 = errorbar(timepoints, squeeze(avisoGRs(i, 1, :)), squeeze(stdevisoGRs(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
%     err2 = errorbar(timepoints, squeeze(avisoGRs(i, 2, :)), squeeze(stdevisoGRs(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
%     err3 = errorbar(timepoints, squeeze(avisoGRs(i, 3, :)), squeeze(stdevisoGRs(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
%     
%     plot(timepoints(privatebenefit(i, 3)), privatebenefit(i, 1), 'x', 'Color', colors(4, :), 'MarkerSize', 10, 'LineWidth', 5)
%     plot(timepoints(privatebenefit(i, 3)), privatebenefit(i, 4), 'x', 'Color', colors(1, :), 'MarkerSize', 10, 'LineWidth', 5)
%     
%     xlabel("Hours");
%     ylabel("GR");
%     xlim([0 24]);
%     ylim([-0.5 1.5]);
%     title(isolatenames(i))
%     legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
%     set(gca,'FontSize',20)
%     axis square
%     hold off;
% end

% % private benefit by crash
% 
% privatebenefit = zeros(nisolates, 5); % GR, OD, timepoint, corresponding noabx GR, calculation
% 
% % first 3 hours, exclude first two timepoint
% [privatebenefit(:, 1), privatebenefit(:, 3)] = min(squeeze(avisoGRs(:, 2, 3:3*6+1)), [], 2);
% for i = 1:nisolates
%     privatebenefit(i, 2) = avisodata(i, 2, privatebenefit(i, 3));
%     privatebenefit(i, 4) = avisoGRs(i, 1, privatebenefit(i, 3));
%     privatebenefit(i, 5) = privatebenefit(i, 1) / privatebenefit(i, 4);
% end
% 
% b = bar(privatebenefit(:, 5));
% b.FaceColor = 'flat';
% b.CData = colors;
% set(gca, 'XTickLabel', isolatenames);
% set(gca, 'FontSize', 20);
% 
% % plot crash
% %GRs
% for i = 1:nisolates   
%     figure(str2double(isolatenames(i))+60000)
%     hold on
%     
%     set1 = plot(timepoints, squeeze(avisoGRs(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
%     set2 = plot(timepoints, squeeze(avisoGRs(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
%     set3 = plot(timepoints, squeeze(avisoGRs(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
%     
%     err1 = errorbar(timepoints, squeeze(avisoGRs(i, 1, :)), squeeze(stdevisoGRs(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
%     err2 = errorbar(timepoints, squeeze(avisoGRs(i, 2, :)), squeeze(stdevisoGRs(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
%     err3 = errorbar(timepoints, squeeze(avisoGRs(i, 3, :)), squeeze(stdevisoGRs(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
%     
%     plot(timepoints(privatebenefit(i, 3)), privatebenefit(i, 1), 'x', 'Color', colors(1, :), 'MarkerSize', 10, 'LineWidth', 5)
%     plot(timepoints(privatebenefit(i, 3)), privatebenefit(i, 4), 'x', 'Color', colors(4, :), 'MarkerSize', 10, 'LineWidth', 5)
%     
%     xlabel("Hours");
%     ylabel("GR");
%     xlim([0 24]);
%     ylim([-0.5 1.5]);
%     title(isolatenames(i))
%     legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
%     set(gca,'FontSize',20)
%     axis square
%     hold off;
% end

% inhibition by abx peak
inhibitioneffect = zeros(nisolates, 5); % GR, OD, timepoint, corresponding inhibited GR, calculation

[inhibitioneffect(:, 1), inhibitioneffect(:, 3)] = max(squeeze(avisoGRs(:, 2, :)), [], 2);
for i = 1:nisolates
    inhibitioneffect(i, 2) = avisodata(i, 2, inhibitioneffect(i, 3));
    inhibitioneffect(i, 4) = avisoGRs(i, 3, inhibitioneffect(i, 3));
    inhibitioneffect(i, 5) = inhibitioneffect(i, 4) / inhibitioneffect(i, 1);
end
% 
% b = bar(inhibitioneffect(:, 5));
% b.FaceColor = 'flat';
% b.CData = colors;
% set(gca, 'XTickLabel', isolatenames);
% set(gca, 'FontSize', 20);
% 
% % plot peak
% %GRs
% for i = 1:nisolates   
%     figure(str2double(isolatenames(i))+70000)
%     hold on
%     
%     set1 = plot(timepoints, squeeze(avisoGRs(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
%     set2 = plot(timepoints, squeeze(avisoGRs(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
%     set3 = plot(timepoints, squeeze(avisoGRs(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
%     
%     err1 = errorbar(timepoints, squeeze(avisoGRs(i, 1, :)), squeeze(stdevisoGRs(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
%     err2 = errorbar(timepoints, squeeze(avisoGRs(i, 2, :)), squeeze(stdevisoGRs(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
%     err3 = errorbar(timepoints, squeeze(avisoGRs(i, 3, :)), squeeze(stdevisoGRs(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
%     
%     plot(timepoints(inhibitioneffect(i, 3)), inhibitioneffect(i, 1), 'x', 'Color', colors(1, :), 'MarkerSize', 10, 'LineWidth', 5)
%     plot(timepoints(inhibitioneffect(i, 3)), inhibitioneffect(i, 4), 'x', 'Color', colors(8, :), 'MarkerSize', 10, 'LineWidth', 5)
%     
%     xlabel("Hours");
%     ylabel("GR");
%     xlim([0 24]);
%     ylim([-0.5 1.5]);
%     title(isolatenames(i))
%     legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
%     set(gca,'FontSize',20)
%     axis square
%     hold off;
% end

% % inhibition by crash
% 
% inhibitioneffect = zeros(nisolates, 5); % GR, OD, timepoint, corresponding abx GR, calculation
% 
% % first five hours, exclude first two timepoint
% [inhibitioneffect(:, 1), inhibitioneffect(:, 3)] = min(squeeze(avisoGRs(:, 3, 3:5*6+1)), [], 2);
% for i = 1:nisolates
%     inhibitioneffect(i, 2) = avisodata(i, 3, inhibitioneffect(i, 3));
%     inhibitioneffect(i, 4) = avisoGRs(i, 2, inhibitioneffect(i, 3));
%     inhibitioneffect(i, 5) = inhibitioneffect(i, 1) / inhibitioneffect(i, 4);
% end
% 
% b = bar(inhibitioneffect(:, 5));
% b.FaceColor = 'flat';
% b.CData = colors;
% set(gca, 'XTickLabel', isolatenames);
% set(gca, 'FontSize', 20);
% 
% % plot crash
% %GRs
% for i = 1:nisolates   
%     figure(str2double(isolatenames(i))+80000)
%     hold on
%     
%     set1 = plot(timepoints, squeeze(avisoGRs(i, 1, :)), 'Color', colors(4, :), 'LineWidth', 3);
%     set2 = plot(timepoints, squeeze(avisoGRs(i, 2, :)), 'Color', colors(1, :), 'LineWidth', 3);
%     set3 = plot(timepoints, squeeze(avisoGRs(i, 3, :)), 'Color', colors(8, :), 'LineWidth', 3);
%     
%     err1 = errorbar(timepoints, squeeze(avisoGRs(i, 1, :)), squeeze(stdevisoGRs(i, 1, :)), 'Color',  colors(4, :) + 0.2*(ones(1, 3) - colors(4, :)), 'CapSize', 0);
%     err2 = errorbar(timepoints, squeeze(avisoGRs(i, 2, :)), squeeze(stdevisoGRs(i, 2, :)), 'Color',  colors(1, :) + 0.2*(ones(1, 3) - colors(1, :)), 'CapSize', 0);
%     err3 = errorbar(timepoints, squeeze(avisoGRs(i, 3, :)), squeeze(stdevisoGRs(i, 3, :)), 'Color',  colors(8, :) + 0.2*(ones(1, 3) - colors(8, :)), 'CapSize', 0);
%     
%     plot(timepoints(inhibitioneffect(i, 3)), inhibitioneffect(i, 1), 'x', 'Color', colors(8, :), 'MarkerSize', 10, 'LineWidth', 5)
%     plot(timepoints(inhibitioneffect(i, 3)), inhibitioneffect(i, 4), 'x', 'Color', colors(1, :), 'MarkerSize', 10, 'LineWidth', 5)
%     
%     xlabel("Hours");
%     ylabel("GR");
%     xlim([0 24]);
%     ylim([-0.5 1.5]);
%     title(isolatenames(i))
%     legend([set1(1), set2(1), set3(1)], 'No antibiotic', 'AMX', 'AMX/CLA', 'Location', 'eastoutside');
%     set(gca,'FontSize',20)
%     axis square
%     hold off;
% end

fig = figure(265928653);
hold on
scatter(privatebenefit(:,5), inhibitioneffect(:, 5), 100, colors, 'filled')
xlabel("private benefit")
ylabel("inhibition effect")
axis square
ylim([0 1.2])
xlim([0, 1.2])

for i = 1:nisolates
    text(privatebenefit(i, 5)+0.02, inhibitioneffect(i, 5), isolatenames(i), 'FontSize', 14)
end
hold off
set(gca, 'FontSize', 20)

saveas(fig, "private benefit vs inhibition", "jpeg");

%% save variables to workspace

save(strcat(isolatenames(1), " ", isolatenames(nisolates)), 'isolatedata', 'isolateGRs', 'isolatenames');
clear
close all

%%

datafile = '4-9-2022 TAZ summary.xlsx';

data = readmatrix(datafile, 'Sheet', 1, 'Range', 'F:I');
info = readmatrix(datafile, 'Sheet', 1, 'Range', 'A:E'); % experiment strain timepoint [a] [i] 

%% split data out by condition
timepoints = unique(info(:, 3));
abxconcentrations = unique(info(~isnan(info(:, 4)), 4));
inhconcentrations = unique(info(~isnan(info(:, 5)), 5));
nstrains = size(unique(info(:, 2)), 1);
nreplicates = size(unique(info(:, 1)), 1); % experiments

% time abx inh strain replicate [fraction OD gfp bfp]
datasplit = nan(length(timepoints), length(abxconcentrations), length(inhconcentrations), nstrains, nreplicates, 4);

for i = 1:size(data, 1)
    datasplit(timepoints==info(i, 3), abxconcentrations==info(i, 4), inhconcentrations==info(i, 5), info(i, 2), info(i, 1), :) = data(i, :);
end
%%

averagedata = squeeze(mean(datasplit, 5, 'omitnan'));
stderrdata = squeeze(std(datasplit, 0, 5, 'omitnan') ./ sqrt(nreplicates));

averagefractions = averagedata(:, :, :, :, 1);
stderrfractions = stderrdata(:, :, :, :, 1);

averageOD = averagedata(:, :, :, :, 2);
stderrOD = stderrdata(:, :, :, :, 2);

averageGFP = averagedata(:, :, :, :, 3);
stderrGFP = stderrdata(:, :, :, :, 3);

averageBFP = averagedata(:, :, :, :, 4);
stderrBFP = stderrdata(:, :, :, :, 4);

GFPperOD = squeeze(datasplit(:, :, :, :, :, 3) ./ datasplit(:, :, :, :, :, 2));
averageGFPperOD = mean(GFPperOD, 5, 'omitnan');
stderrGFPperOD = std(GFPperOD, 0, 5, 'omitnan') ./ sqrt(nreplicates);
%% plot

% time abx inh strain replicate [fraction OD GFP BFP]

strainindex = 1;
strainnames = ["Low Copy Bla", "Low Copy BlaM", "High Copy Bla", "High Copy BlaM"];

%% average with stderr

for i = 1:nstrains
    for j = 2
    %for j = 1:length(abxconcentrations)
        f = figure(10*i + j);
        hold on
        errorbar(inhconcentrations, squeeze(averagefractions(2, j, :, i)), squeeze(stderrfractions(2, j, :, i)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
        %yline(averagefractions(1, 1, 1, i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        axis square
        ylim([0 1])
        xlim([0 max(inhconcentrations)])
        %xlabel("TAZ (\mug/mL)")
        %ylabel("Resistant Fraction")
        
        %title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX"))
        set(gca, 'FontSize', 20)
        
        %%saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX Fraction Averages"), 'jpeg')
    end
end

close all
%% plot Bla and BlaM together

% % low copy
% figure(100)
% hold on
% errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 1)), squeeze(stderrfractions(2, 2, :, 1)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5)
% errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 2)), squeeze(stderrfractions(2, 2, :, 2)), 'o-', 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'MarkerSize', 10, 'LineWidth', 1.5)
% axis square
% ylim([0 1])
% xlim([0 max(inhconcentrations)])
% %xlabel("TAZ (\mug/mL)")
% %ylabel("Resistant Fraction")
% %legend('Periplasmic', 'Cytoplasmic', 'Location', 'northeast', 'FontSize', 24)
% %legend('boxoff')
% set(gca, 'FontSize', 30)
% 
% %high copy
% figure(101)
% hold on
% errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 3)), squeeze(stderrfractions(2, 2, :, 3)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5)
% errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 4)), squeeze(stderrfractions(2, 2, :, 4)), 'o-', 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'MarkerSize', 10, 'LineWidth', 1.5)
% axis square
% ylim([0 1])
% xlim([0 max(inhconcentrations)])
% %xlabel("TAZ (\mug/mL)")
% %ylabel("Resistant Fraction")
% %legend('Bla', 'BlaM', 'Location', 'eastoutside')
% set(gca, 'FontSize', 30)
% 
% %close all

%% plot with replicates

% low copy
figure(200)
hold on
% scatter replicates
scatter(inhconcentrations, squeeze(datasplit(2, 2, :, 1, :, 1)), 100, [0.7 0.7 0.7], 'filled');
scatter(inhconcentrations, squeeze(datasplit(2, 2, :, 2, :, 1)), 100, [0.2 0.2 0.2]);
% plot averages
errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 1)), squeeze(stderrfractions(2, 2, :, 1)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5)
errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 2)), squeeze(stderrfractions(2, 2, :, 2)), 'o-', 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'MarkerSize', 10, 'LineWidth', 1.5)
axis square
ylim([0 1])
xlim([0 max(inhconcentrations)])
%xlabel("TAZ (\mug/mL)")
%ylabel("Resistant Fraction")
%legend('Periplasmic', 'Cytoplasmic', 'Location', 'northeast', 'FontSize', 24)
%legend('boxoff')
set(gca, 'FontSize', 30)

%high copy
figure(201)
hold on
% scatter replicates
scatter(inhconcentrations, squeeze(datasplit(2, 2, :, 3, :, 1)), 100, [0.7 0.7 0.7], 'filled');
scatter(inhconcentrations, squeeze(datasplit(2, 2, :, 4, :, 1)), 100, [0.2 0.2 0.2]);
% plot averages
errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 3)), squeeze(stderrfractions(2, 2, :, 3)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5)
errorbar(inhconcentrations, squeeze(averagefractions(2, 2, :, 4)), squeeze(stderrfractions(2, 2, :, 4)), 'o-', 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'MarkerSize', 10, 'LineWidth', 1.5)
axis square
ylim([0 1])
xlim([0 max(inhconcentrations)])
%xlabel("TAZ (\mug/mL)")
%ylabel("Resistant Fraction")
%legend('Bla', 'BlaM', 'Location', 'eastoutside')
set(gca, 'FontSize', 30)

%close all


%% OD 
% 
% % replicates
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j);
%         for k = 1:nreplicates
%             plot(inhconcentrations, squeeze(datasplit(2, j, :, i, k, 2)), 'o-', 'MarkerFaceColor', 'auto') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%             hold on
%         end
%         axis square
%         ylim([0 1.5])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("OD600")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX OD Replicates"), 'jpeg')
%     end
% end
% 
% close all

% average with stderr
% 
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j + 100);
%         hold on
%         errorbar(inhconcentrations, squeeze(averageOD(2, j, :, i)), squeeze(stderrOD(2, j, :, i)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%         % yline(averageOD(1, 1, 1, i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
%         axis square
%         ylim([0 1.5])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("OD600")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX OD"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX OD Averages"), 'jpeg')
%     end
% end
% 
% close all
% 

%% GFP
% 
% % replicates
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j);
%         for k = 1:nreplicates
%             plot(inhconcentrations, squeeze(datasplit(2, j, :, i, k, 3)), 'o-', 'MarkerFaceColor', 'auto') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%             hold on
%         end
%         axis square
%         ylim([0 5000])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("GFP")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX GFP"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX GFP Replicates"), 'jpeg')
%     end
% end
% 
% close all

% average with stderr
% 
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j + 100);
%         hold on
%         errorbar(inhconcentrations, squeeze(averageGFP(2, j, :, i)), squeeze(stderrGFP(2, j, :, i)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%         % yline(averageOD(1, 1, 1, i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
%         axis square
%         ylim([0 5000])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("GFP")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX GFP"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX GFP Averages"), 'jpeg')
%     end
% end
% 
% close all

%% BFP
% 
% % replicates
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j);
%         for k = 1:nreplicates
%             plot(inhconcentrations, squeeze(datasplit(2, j, :, i, k, 4)), 'o-', 'MarkerFaceColor', 'auto') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%             hold on
%         end
%         axis square
%         ylim([0 500])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("BFP")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX BFP"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX BFP Replicates"), 'jpeg')
%     end
% end
% 
% close all

% average with stderr
% 
% for i = 1:nstrains
%     for j = 1:length(abxconcentrations)
%         f = figure(10*i + j + 100);
%         hold on
%         errorbar(inhconcentrations, squeeze(averageBFP(2, j, :, i)), squeeze(stderrBFP(2, j, :, i)), 'o-', 'Color', 'k', 'MarkerFaceColor', 'k') %'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2])
%         % yline(averageOD(1, 1, 1, i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
%         axis square
%         ylim([0 500])
%         xlim([0 max(inhconcentrations)])
%         xlabel("TAZ (\mug/mL)")
%         ylabel("BFP")
%         
%         title(strcat(strainnames(i), " + ", num2str(abxconcentrations(j)), " \mug/mL AMX BFP"))
%         set(gca, 'FontSize', 20)
%         
%         %saveas(f, strcat(strainnames(i), " ", num2str(abxconcentrations(j)), " AMX BFP Averages"), 'jpeg')
%     end
% end
% 
% close all
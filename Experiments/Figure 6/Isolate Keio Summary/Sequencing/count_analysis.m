clear
close all

%%
% read in raw data from excel
datafile = "merged_counts_20240325";
optsinfo = detectImportOptions(datafile, 'Sheet', 1, 'Range', 'A1:H193');
vartypes = optsinfo.VariableTypes;
vartypes{1} = 'string';
vartypes{8} = 'string';
optsinfo.VariableTypes = vartypes;
info = readtable(datafile,optsinfo);

data = readmatrix(datafile, 'Sheet', 1, 'Range', 'I2:AO193');

[nwells, nkeiostrains] = size(data);
nreplicates = max(info{:, 7});
rstrainnames = unique(info{:, 1});
timepoints = unique(info{:, 3});
treatments = unique(info{:, 4:6}, 'rows');
ntreatments = length(treatments);
keiostrainnames = strcat("K", string(readmatrix(datafile, 'Range', 'I1:AO1')));

%% sum totals and normalize

totals = sum(data, 2);
normalized = data ./ totals;

%% 

colors = cmocean('diff', nkeiostrains);
%colors = crameri("oslo", nstrains);

% optional to increase visual changes
% permutevector = randperm(nstrains);
% colors = colors(permutevector, :);

% cmocean thermal
% fcrameri romaO

figure(1082357)
barcodes = info{:, 8};
b = bar(barcodes, normalized, 'stacked', 'FaceColor', 'flat');
for i = 1:nkeiostrains
    b(i).CData = colors(i, :);
end
leg = legend(keiostrainnames, "Location", 'eastoutside');
set(leg, 'FontSize', 12)
box(leg, "off")
ylim([0 1])

%% split into replicates

conditiondata = zeros(6, 9, 2, 5, 3, nkeiostrains); % Keio set, R strain, timepoint, 5 abx conditions, replicates, strains
normalizedconditiondata = conditiondata;

for i = 1:nwells
    [~, ~, treatmentid] = intersect(info{i, 4:6}, treatments, 'rows');
    conditiondata(info{i, 2},  info{i, 1} == rstrainnames, info{i, 3} == timepoints, treatmentid, info{i, 7}, :) = data(i, :);
    normalizedconditiondata(info{i, 2},  info{i, 1} == rstrainnames, info{i, 3} == timepoints, treatmentid, info{i, 7}, :) = normalized(i, :);
end

%% plot by condition

% for i = 1:length(rstrainnames)-1
%     switch i
%         case {1, 2}
%             keio = 1;
%         case {3, 4}
%             keio = 2;
%         case {5, 6}
%             keio = 3;
%         case {7, 8}
%             keio =4;
%     end
%     for j = 1:ntreatments
%         f = figure(100*i + j);
%         f.Position = [100 100 180 500];
%         repb = bar(squeeze(normalizedconditiondata(keio, i, 2, j, :, :)), 0.85, 'stacked', 'FaceColor', 'flat');
%         for m = 1:nkeiostrains
%             repb(m).CData = colors(m, :);
%         end
%         ylim([0 1])
%     end
% end

%% average normalized replicates (average fraction)
avgnormalizedconditiondata = squeeze(mean(normalizedconditiondata, 5)); % Keio set, R strain, timepoint, ABX condition, strains

%% plot avgs
% with resistant strains
for i = 1:length(rstrainnames)-1
    switch i
        case {1, 2}
            keio = 1;
        case {3, 4}
            keio = 2;
        case {5, 6}
            keio = 3;
        case {7, 8}
            keio =4;
    end
    f = figure(100+i);
    hold on
    f.Position = [100 100 220 500];
    %repb = bar(squeeze(avgnormalizedconditiondata(keio, i, 2, :, :)), 0.9, 'stacked', 'FaceColor', 'flat');
    repb = bar([0, 2, 4], squeeze(avgnormalizedconditiondata(keio, i, 2, 3:5, :)), 0.85, 'stacked', 'FaceColor', 'flat');
    for m = 1:nkeiostrains
        repb(m).CData = colors(m, :);
    end
    ylim([0 1])
    % xlabel("[CLA] \mug/mL")
    % title(rstrainnames(i))
    set(gca, 'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1])
    set(gca, 'XTick', [])
    box off
    axis off
    bl = repb.BaseLine;
    bl.Visible = 'off';
    saveas(f, strcat("Figure 6B ", rstrainnames(i)), 'svg');
end

% keio
for i = 5
    f = figure(1000+i);
    hold on
    f.Position = [100 100 220 500];
    %repb = bar(squeeze(avgnormalizedconditiondata(i, 9, 2, :, :)), 0.9, 'stacked', 'FaceColor', 'flat');
    repb = bar([0, 2, 4], squeeze(avgnormalizedconditiondata(i, 9, 2, 3:5, :)), 0.85, 'stacked', 'FaceColor', 'flat');
    for m = 1:nkeiostrains
        repb(m).CData = colors(m, :);
    end
    ylim([0 1])
    % xlabel("[CLA] \mug/mL")
    % title("Pure Keio")
    set(gca, 'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1])
    set(gca, 'XTick', [])
    box off
    axis off
    bl = repb.BaseLine;
    bl.Visible = 'off';
    saveas(f, "Figure 6B Keio", 'svg');
end

%% plot initial timepoints
% with resistant strains
for i = 1:length(rstrainnames)-1
    switch i
        case {1, 2}
            keio = 1;
        case {3, 4}
            keio = 2;
        case {5, 6}
            keio = 3;
        case {7, 8}
            keio =4;
    end
    f = figure(200+i);
    f.Position = [100 100 220 500];
    %repb = bar(squeeze(avgnormalizedconditiondata(keio, i, 2, :, :)), 0.9, 'stacked', 'FaceColor', 'flat');
    repb = bar([0, 2, 4], squeeze(avgnormalizedconditiondata(keio, i, 1, 1:3, :)), 0.85, 'stacked', 'FaceColor', 'flat');
    for m = 1:nkeiostrains
        repb(m).CData = colors(m, :);
    end
    ylim([0 1])
    %xlabel("[CLA] \mug/mL")
    % title(rstrainnames(i))
    set(gca, 'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1])
    set(gca, 'XTick', [])
    box off
    axis off
    bl = repb.BaseLine;
    bl.Visible = 'off';
    saveas(f, strcat("Figure 6B ", rstrainnames(i), " initial"), 'svg');
end

% keios
for i = 5
    f = figure(2000+i);
    f.Position = [100 100 220 500];
    %repb = bar(squeeze(avgnormalizedconditiondata(i, 9, 2, :, :)), 0.9, 'stacked', 'FaceColor', 'flat');
    repb = bar([0, 2, 4], squeeze(avgnormalizedconditiondata(i, 9, 1, 1:3, :)), 0.85, 'stacked', 'FaceColor', 'flat');
    for m = 1:nkeiostrains
        repb(m).CData = colors(m, :);
    end
    ylim([0 1])
    % xlabel("[CLA] \mug/mL")
    % title("Pure Keio")
    set(gca, 'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1])
    set(gca, 'XTick', [])
    box off
    axis off
    bl = repb.BaseLine;
    bl.Visible = 'off';
    saveas(f, strcat("Figure 6B Keio initial"), 'svg');
end
%% for each strain, get enrichment vs initial

% %averagenormalized dimensions: Keio set, R strain, timepoint, 5 abx conditions, strains
% enrichmentfactors = zeros(6, 9, 5, nkeiostrains);
% foldenrichment = enrichmentfactors;
% 
% for i = 1:6
%     for j = 1:length(rstrainnames)
%         for k = 1:ntreatments
%             [~, ~, zerocheck] = intersect(squeeze(avgnormalizedconditiondata(i, j, 2, k, :))', zeros(1, nkeiostrains), 'rows');
%             if isempty(zerocheck)
%                 final = squeeze(avgnormalizedconditiondata(i, j, 2, k, :));
%                 initial = squeeze(avgnormalizedconditiondata(i, j, 1, 1, :));
%                 enrichmentfactors(i, j, k, :) =  final ./ initial;
%                 foldenrichment(i, j, k, :) = (final - initial) ./ initial;
%             end
%         end
%     end
% end

%% plot 

% colors = crameri('-devon', ntreatments);
% for i = 1:nkeiostrains
%     % should be 8 strains and 2 keio sets with 5 conditions each
%     longEFs = zeros(10, 5);
%     for j = 1:10
%         rstrain = j;
%         switch j
%             case {1, 2}
%                 keio = 1;
%             case {3, 4}
%                 keio = 2;
%             case {5, 6}
%                 keio = 3;
%             case {7, 8}
%                 keio =4;
%             case 9
%                 keio = 5;
%                 rstrain = 9;
%             case 10
%                 keio = 6;
%                 rstrain = 9;
%         end
%         longEFs(j, :) = squeeze(enrichmentfactors(keio, rstrain, :, i));
%         longfoldEFs(j, :) = squeeze(foldenrichment(keio, rstrain, :, i));
%     end
% 
%     f = figure(10000+i);
%     hold on
%     kstrainb = bar(longfoldEFs, 'FaceColor', 'flat');
%     for j = 1:ntreatments
%         kstrainb(j).CData = colors(j, :)
%     end
%     ylim([-1 5])
%     ylabel("Fold Enrichment")
%     set(gca, 'XTick', 1:10)
%     set(gca, 'XTickLabels', cat(1, rstrainnames(1:8), ["Pure Keio", "Pure Keio"]'))
%     set(gca, 'FontSize', 15)
%     axis square
% end

%% plot in betamin order

% order = [4 1 7 5 8 6 2 3 9 10];
% 
% colors = crameri('-devon', ntreatments);
% 
% for i = 1%:nkeiostrains
%     % should be 8 strains and 2 keio sets with 5 conditions each
%     longEFs = zeros(10, 5);
%     for j = 1:10
%         rstrain = order(j);
%         switch rstrain
%             case {1, 2}
%                 keio = 1;
%             case {3, 4}
%                 keio = 2;
%             case {5, 6}
%                 keio = 3;
%             case {7, 8}
%                 keio =4;
%             case 9
%                 keio = 5;
%                 rstrain = 9;
%             case 10
%                 keio = 6;
%                 rstrain = 9;
%         end
%         longEFs(j, :) = squeeze(enrichmentfactors(keio, rstrain, :, i));
%         longfoldEFs(j, :) = squeeze(foldenrichment(keio, rstrain, :, i));
%     end
% 
%     f = figure(20000+i);
%     hold on
%     orderedkstrainb = bar(longfoldEFs, 'FaceColor', 'flat');
%     %orderedkstrainb = bar(longfoldEFs(:, 3:5), 'FaceColor', 'flat');
%     for j = 1:ntreatments
%     %for j = 1:3
%         orderedkstrainb(j).CData = colors(j, :);
%         %orderedkstrainb(j).CData = colors(j+2, :);
%     end
%     ylim([-1 5])
%     ylabel("Fold Enrichment")
%     set(gca, 'XTick', 1:10)
%     set(gca, 'XTickLabels', cat(1, rstrainnames(order(1:8)), ["Pure Keio", "Pure Keio"]'))
%     set(gca, 'FontSize', 15)
%     axis square
%     title(keiostrainnames(i))
%     l = legend(["untreated", "cm only", "CLA 0", "CLA 2", "CLA 4"], 'Location', 'eastoutside');
%     % l = legend(string([0, 2, 4]), 'Location', 'eastoutside');
%     % title(l, "[CLA] (\mug/mL)")
%     box(l, "off")
%     %saveas(f, keiostrainnames(i), 'jpg')
% end

%% betamin order tiled layout

% order = [4 1 7 5 8 6 2 3 9 10];
% 
% colors = crameri('-devon', ntreatments);
% 
% f = figure(298572);
% tiledlayout(4, 9, 'TileSpacing', 'tight')
% counter = 1;
% for i = 1:nkeiostrains
%     % should be 8 strains and 2 keio sets with 5 conditions each
%     longEFs = zeros(10, 5);
%     for j = 1:10
%         rstrain = order(j);
%         switch rstrain
%             case {1, 2}
%                 keio = 1;
%             case {3, 4}
%                 keio = 2;
%             case {5, 6}
%                 keio = 3;
%             case {7, 8}
%                 keio =4;
%             case 9
%                 keio = 5;
%                 rstrain = 9;
%             case 10
%                 keio = 6;
%                 rstrain = 9;
%         end
%         longEFs(j, :) = squeeze(enrichmentfactors(keio, rstrain, :, i));
%         longfoldEFs(j, :) = squeeze(foldenrichment(keio, rstrain, :, i));
%     end
% 
%     nexttile
%     orderedkstrainb = bar(longfoldEFs, 'FaceColor', 'flat');
%     %orderedkstrainb = bar(longfoldEFs(:, 3:5), 'FaceColor', 'flat');
%     for j = 1:ntreatments
%     %for j = 1:3
%         orderedkstrainb(j).CData = colors(j, :);
%         %orderedkstrainb(j).CData = colors(j+2, :);
%     end
%     ylim([-1 5])
%     set(gca, 'XTick', [])
%     if counter == (9*3+1)
%         set(gca, 'XTick', 1:10)
%         set(gca, 'XTickLabels', cat(1, rstrainnames(order(1:8)), ["Pure Keio", "Pure Keio"]'))
%         ylabel("Fold Enrichment")
%     end
%     set(gca, 'FontSize', 12)
%     axis square
%     title(keiostrainnames(i))
%     counter = counter + 1;
% end
% %close all
%% batch_data_analysis.m
% combine isolate data into one massive matrix, analyze. When adding new
% data, enter Y at the prompt and then type the name of the .mat file.
% bug: the rest of the code won't run properly when you're adding new data.
% Just ignore, clear, and rerun saying N. 

% isolate x treatment x replicates (12, 4 tech reps for each bio rep in order) x timepoints

clear
close all

%% add new data to old combined data, one at a time

load("combineddata")
load("anonisolatenames") 

%% set constants

% set or get constants
nwells = 384;
height = 16;
width = 24;
ntimepoints = 145;
nisolates = length(combinedisolatenames);
ntreatments = 3;
ntechreps = 4;
nbioreps = 3;
ntotalreps = ntechreps * nbioreps;
timepoints = linspace(0, 24, ntimepoints);

%% Figure 5A

colors = [215 48 39;
        244 109 67;
        253 174 97;
        254 224 144;
        224 243 248;
        171 217 233;
        116 173 209;
        69 117 180];
colors = colors ./ 255;

fig = figure(34234324);
counter = 1;

shift=0;
indices = [56 87 104 172 202 204 216 247];

plotteddata = zeros(nreplicates*nconditions*8, 145+3);

for i = indices(:, [1 7 4 5 6 8 2 3])
    
    % combinedisolate data: isolate, condition, replicate, timepoints
    %figure(i)
    subplot(2, 4, counter);
    hold on
    
    set1 = plot(timepoints, squeeze(combinedisolatedata(i+shift, 1, :, :)), 'Color', colors(4, :), 'LineWidth', 3);
    set2 = plot(timepoints, squeeze(combinedisolatedata(i+shift, 2, :, :)), 'Color', colors(1, :), 'LineWidth', 3);
    set3 = plot(timepoints, squeeze(combinedisolatedata(i+shift, 3, :, :)), 'Color', colors(8, :), 'LineWidth', 3);

    start = (counter-1)*nreplicates*nconditions+1;
    plotteddata(start:start+nreplicates*nconditions-1, 1) = repelem(combinedisolatenames(i), 36);
    plotteddata(start:start+nreplicates*nconditions-1, 2) = repelem(1:3, nreplicates);
    plotteddata(start:start+nreplicates*nconditions-1, 3) = repmat(1:12, 1, 3);
    plotteddata(start:start+nreplicates-1, 4:end) = squeeze(combinedisolatedata(i+shift, 1, :, :));
    plotteddata(start+nreplicates:start+2*nreplicates-1, 4:end) = squeeze(combinedisolatedata(i+shift, 2, :, :));
    plotteddata(start+2*nreplicates:start+3*nreplicates-1, 4:end) = squeeze(combinedisolatedata(i+shift, 3, :, :));

    if counter == 5
        xlabel("Hours");
        ylabel("OD");
    else
        set(gca, 'YTickLabel', [], 'XTickLabel', []);
    end
    xlim([0 24]);
    ylim([0 2]);
    title(combinedisolatenames(i))
    legend([set1(1), set2(1), set3(1)], "No antibiotic", "50 \mug/mL AMX", "50 \mug/mL AMX +" + newline + "25 \mug/mL CLA", 'Location', 'eastoutside'); 
    legend('boxoff')
    set(gca, 'XTick', linspace(0, 24, 3))
    set(gca,'FontSize',20)
    axis square
    hold off;
    counter = counter + 1;


end

%% get averages
avisoGRs = squeeze(mean(combinedisolateGRs, 3)); % isolate x treatment (no, amx, amxcla) x timepoint
avisodata = squeeze(mean(combinedisolatedata, 3)); % isolate x treatment (no, amx, amxcla) x timepoint

stdevisoGRs = squeeze(std(combinedisolateGRs, 0, 3)) / sqrt(ntotalreps); % 0 is default weight
stdevisodata = squeeze(std(combinedisolatedata, 0, 3)) / sqrt(ntotalreps);

%% analyze private benefit stuff

% private benefit by crash (resistance)

privatebenefit = zeros(nisolates, 5); % GR, OD, index, corresponding noabx GR, calculated value

% first 3 hours, exclude first two timepoint

skip = 2;
[privatebenefit(:, 1), privatebenefit(:, 3)] = min(squeeze(avisoGRs(:, 2, 1+skip:3*6+1)), [], 2);
for i = 1:nisolates
    privatebenefit(i, 2) = avisodata(i, 2, privatebenefit(i, 3)+skip);
    privatebenefit(i, 4) = avisoGRs(i, 1, privatebenefit(i, 3)+skip);
    privatebenefit(i, 5) = privatebenefit(i, 1) / privatebenefit(i, 4);
end

% inhibition by crash

inhibitioneffect = zeros(nisolates, 5); % GR, OD, index, corresponding abx GR, calculated value

% first five hours, exclude first two timepoint
skip = 2;

[inhibitioneffect(:, 1), inhibitioneffect(:, 3)] = min(squeeze(avisoGRs(:, 3, 1+skip:5*6+1)), [], 2);
for i = 1:nisolates
    inhibitioneffect(i, 2) = avisodata(i, 3, inhibitioneffect(i, 3)+skip);
    inhibitioneffect(i, 4) = avisoGRs(i, 1, inhibitioneffect(i, 3)+skip);
    inhibitioneffect(i, 5) = inhibitioneffect(i, 1) / inhibitioneffect(i, 4);
end

% recovery time vs no antibiotic
RT = zeros(nisolates, 8); % max noabx OD, thresholdOD, noabx index, noabx OD, amx index, amx OD, amxcla index, amxcla OD

RT(:, 1) = max(squeeze(avisodata(:, 1, :)), [], 2);
RT(:, 2) = 0.5 * RT(:, 1);

for i = 1:nisolates
    for j = 1:ntreatments
        thisisodata = squeeze(avisodata(i, j, :));
        idx = find(thisisodata >= RT(i, 2), 1);
        if ~isempty(idx) 
            RT(i, 2*j+1) = idx;
            RT(i, 2*j+2) = thisisodata(idx);
        else
            RT(i, 2*j+1) = Inf;
            RT(i, 2*j+2) = NaN;
        end
    end
end

resilience = zeros(nisolates, 2); % amx resilience, amxcla resilience

resilience(:, 1) = RT(:, 3) ./ RT(:, 5);
resilience(:, 2) = RT(:, 3) ./ RT(:, 7);
%% Figure 5B
fig = figure(1);
hold on
% color by amx resilience
map = cmocean('amp');
colors = zeros(nisolates, 3);
for i = 1:nisolates
    idx = round((resilience(i, 1) ./ max(resilience(:, 1))).*256);
    if (idx == 0) 
        idx = 1;
    end
    colors(i, :) = map(idx, :);
end

scatter(privatebenefit(:,5), inhibitioneffect(:, 5), 100, colors, 'filled')
colormap(map)
colorbar
% xlabel("resistance")
% ylabel("inhibition effect by crash (reversed)")
axis square
ylim([-0.8 1.5])
xlim([-0.8 1.5])

hold off
set(gca, 'FontSize', 20)
xlabel("AMX resistance")
ylabel("AMX + CLA resistance")

% saveas(fig, "private benefit vs inhibition", "jpeg");

%% beta

%sensitiveindices = find(resilience(:, 1) == 0);
[minABXGRs(:, 1), minABXGRs(:, 2)] = min(squeeze(avisoGRs(:, 2, :)), [], 2);
[peak0GRs(:, 1), peak0GRs(:, 2)] = max(squeeze(avisoGRs(:, 1, :)), [], 2);

% beta * gamma * g = rho(noAB) - rho(maxAB)
[minAIGRs(:, 1), minAIGRs(:, 2)] = min(squeeze(avisoGRs(:, 3, :)), [], 2);

betas = peak0GRs(:, 1) - minAIGRs(:, 1);
figure(18735927)
histogram(betas)
axis square
set(gca, 'FontSize', 20)
ylabel("Number of isolates")
xlabel("\beta")

%% long data
nstrains = 311;
nconditions = 3;
nreplicates = 12;
ntimepoints = 145; 
timepoints = (0:ntimepoints-1)./6;
amxconcs = [0 50 50];
claconcs = [0 0 25];

longdata = nan(nstrains * nconditions * nreplicates * ntimepoints, 5);
%amx concentration, cla concentration, replicate, timepoint, value

for i = 1:nstrains
    for j = 1:nconditions
        for k = 1:nreplicates
            for m = 1:ntimepoints
                index = (i-1)*(nconditions*nreplicates*ntimepoints)+ (j-1)*(nreplicates*ntimepoints) + (k-1)*ntimepoints + m;
                longdata(index, :) = [amxconcs(j), claconcs(j), k, timepoints(m), squeeze(combinedisolatedata(i, j, k, m))];
            end
        end
    end
end

allnames = repelem(combinedisolatenames, nconditions*nreplicates*ntimepoints);

alldata = [table(allnames', 'VariableNames', ["Strain"]) array2table(longdata, 'VariableNames', ["[AMX]", "[CLA]", "Replicate", "Hours", "OD600"])];
%writetable(alldata, "allisolates.csv")
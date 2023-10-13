%% mantis_generate_7282018.m

% designs conditions and scrambles them for plate layout

%% Plate Information

% plate dimensions
height = 32;
width = 48;
nblanks = 128;
ncontrols = 128;
nwells = height * width - (nblanks + ncontrols);
nreplicates = 5;
nconditions = nwells / nreplicates;
% conditions = 256 = 16 * 16
% logspace of 15 points plus 0

%% Reagent Details

antibiotic = [0 logspace(log10(1), log10(500), 15)]; % ug/mL 
inhibitor = [0 logspace(log10(0.05), log10(5), 15)]; % ug/mL

% stock concentrations (ug/uL)
stock_conc_antibiotic_high = 2; % can choose these. Regular stock is 100 mg/mL. Dilute 20 in 1000.
stock_conc_inhibitor_high = 0.02; % can choose these. Regular stock is 2 mg/mL. Dilute 10 in 1000.
stock_conc_antibiotic_low = stock_conc_antibiotic_high ./ 25; % Dilute 40 in 1000.
stock_conc_inhibitor_low = stock_conc_inhibitor_high ./ 25; % Dilute 40 in 1000.

% ug of compound needed in 10 uL culture
desired_ug_antibiotic = antibiotic ./ 100;
desired_ug_inhibitor = inhibitor ./ 100;

% uL of stock to add to 10 uL culture
stock_toadd_antibiotic_low = desired_ug_antibiotic ./ stock_conc_antibiotic_low;
stock_toadd_inhibitor_low = desired_ug_inhibitor ./ stock_conc_inhibitor_low;
stock_toadd_antibiotic_high = desired_ug_antibiotic ./ stock_conc_antibiotic_high;
stock_toadd_inhibitor_high = desired_ug_inhibitor ./ stock_conc_inhibitor_high;

max_toadd = 2.4;

% divide into high and low concentrations
overflow_indices_antibiotic = stock_toadd_antibiotic_low > max_toadd;
stock_toadd_antibiotic_low(overflow_indices_antibiotic) = 0;
stock_toadd_antibiotic_high(~overflow_indices_antibiotic) = 0;
overflow_indices_inhibitor = stock_toadd_inhibitor_low > max_toadd;
stock_toadd_inhibitor_low(overflow_indices_inhibitor) = 0;
stock_toadd_inhibitor_high(~overflow_indices_inhibitor) = 0;

% round off for Mantis
toadd_antibiotic_low = round(round(stock_toadd_antibiotic_low, 2), 1);
toadd_antibiotic_high = round(round(stock_toadd_antibiotic_high, 2), 1);
toadd_inhibitor_low = round(round(stock_toadd_inhibitor_low, 2), 1);
toadd_inhibitor_high = round(round(stock_toadd_inhibitor_high, 2), 1);

%% create matrices

conditionmap = [(1:nconditions).', zeros(nconditions, 4)];
for i = 1:length(antibiotic)
    for j = 1:length(inhibitor)
        index = (i - 1)*length(inhibitor) + j;
        conditionmap(index, 2) = toadd_antibiotic_low(i);
        conditionmap(index, 3) = toadd_antibiotic_high(i);
        conditionmap(index, 4) = toadd_inhibitor_low(j);
        conditionmap(index, 5) = toadd_inhibitor_high(j);
    end
end

% make and shuffle condition matrix
condition_replicates = repmat(1:nconditions, 1, nreplicates);
% add blanks and controls
conditions_withcontrols = cat(2, condition_replicates, zeros(1, ncontrols));
conditions_withblanks = cat(2, conditions_withcontrols, (-1 * ones(1, nblanks)));
conditions_shuffled = conditions_withblanks(randperm(size(conditions_withblanks, 2)));
conditions = reshape(conditions_shuffled, height, width);

%% make maps for mantis
culture = double(conditions >= 0); % all wells except blanks get cells, 1:1000 dilution

antibiotic_low = zeros(height, width);
antibiotic_high = zeros(height, width);
inhibitor_low = zeros(height, width);
inhibitor_high = zeros(height, width);
for i = 1:height
    for j = 1:width
        if conditions(i, j) > 0 % blanks and controls have neither antibiotic nor inhibitor
            antibiotic_low(i, j) = conditionmap(conditions(i, j), 2);
            antibiotic_high(i, j) = conditionmap(conditions(i, j), 3);
            inhibitor_low(i, j) = conditionmap(conditions(i, j), 4);
            inhibitor_high(i, j) = conditionmap(conditions(i, j), 5);
        end
    end
end

medium = zeros(height, width);
medium = 10 - (culture + antibiotic_low + antibiotic_high + inhibitor_low + inhibitor_high);
medium_low = mod(medium, 1);
medium_high = medium - medium_low;

%% Backcalculate actual concentrations (in ug/mL)

conditions_concentrations = [(1:nconditions).', zeros(nconditions, 2)];
for i = 1:nconditions
        conditions_concentrations(i, 2) = 100 * (conditionmap(i, 2) * stock_conc_antibiotic_low ...
            + conditionmap(i, 3) * stock_conc_antibiotic_high);
        conditions_concentrations(i, 3) = 100 * (conditionmap(i, 4) * stock_conc_inhibitor_low ...
            + conditionmap(i, 5) * stock_conc_inhibitor_high);
end

%%

% total volumes

volume_antibiotic_low = sum(sum(antibiotic_low));
volume_antibiotic_high = sum(sum(antibiotic_high));
volume_inhibitor_low = sum(sum(inhibitor_low));
volume_inhibitor_high = sum(sum(inhibitor_high));
volume_medium_low = sum(sum(medium_low));
volume_medium_high = sum(sum(medium_high));
volume_culture = sum(sum(culture));

%% Vikas Pejaver
% Icahn School of Medicine at Mount Sinai
% 2023-2024

%% Wrapper script to add thresholds for 3 points as per ACMG/AMP SVC4
% Also cleans up some of the older files in which multiple tool outputs
% were combined into a single file by separating them for each tool
% Input files: MAT file output by main.m in 
%              'local_posterior_probability'. Files located in
%              'results -> bootstrapped'
% Output file: Same MAT file with updated posterior probability and score
%              threshold information

%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

%% Constants and defaults
in_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/alphamissense/esm1b_10k_3prct_100pts.mat'; %'/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts.mat';
out_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/alphamissense/w3pt_esm1b_10k_3prct_100pts.mat'; %'/Users/vikaspejaver/Desktop/ClinGen-SVI/results/alphamissense/w3pt_vest4_10k_3prct_100pts.mat';
orig_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/data/alphamissense/ESM1b_score_PLP_BLB_predictions.txt'; %'/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/VEST4_score_PLP_BLB_predictions.txt'; % original PLP_BLB file used to compute thresholds; if empty, thresholds will not be recomputed (HAD TO INTRODUCE THIS HACK AS WE DIDN'T STORE TOOL_SPECIFIC THRESHOLD RANGES BEFORE)
method_name = 'ESM1b'; %'VEST4';
i = 1; %11; % tool index in input file
to_drop = [2:4]; % posterior thresholds to drop


%% Load MAT file and initialize/update relevant variables for output
currvars = load(in_file);
finalvars = currvars;
finalvars.DiscountedThresholdP = {};
finalvars.DiscountedThresholdB = {};
finalvars.ThresholdP = {};
finalvars.ThresholdB = {};
finalvars.pthresh = {};
finalvars.bthresh = {};
finalvars.posteriors_p = {};
finalvars.posteriors_b = {};


%% If unique thresholds need to be recomputed, do so
if ~strcmp(orig_file, '')
    % read in PLP/BLB variants
    D = load(orig_file);
    x = D(:, 1);
    y = D(:, 2);

    % read in gnomAD variants
    file = strrep(orig_file, 'BLB', 'U');
    D = load(file);
    g = D(D(:, 2) == 0, 1);

    % negate the scores for methods where high prediction means benign
    if strcmp(method_name, 'SIFT') == 1 || strcmp(method_name, 'FATHMM') == 1 || strcmp(method_name, 'ESM1b') == 1
        x = -x;
        g = -g;
    end

    % must remove one gnomAD variant with score > 1, as per Panos
    %if strcmp(method_name, 'EA1.0') == 1
    %    g = g(g <= 1);
    %end

    % the worth of a negative example to satisfy the prior alpha
    w = (1 - finalvars.alpha) * sum(y == 1) / (sum(y == 0) * finalvars.alpha);

    % thresholds for posterior for pathogenicity, reverse sorted
    thrs = flip(unique([x; g; floor(min([x; g])); ceil(max([x; g]))]));
end


%% Thresholds for posteriors, p and b, for alpha given above
% 8 pts = VS, 4 pts = ST, 2 pts = MO, 1 pts = SU (but also computes for
% every integer point value in between
Post_p = [];
Post_b = [];
for j = 8 : -1 : 1
    Post_p = [Post_p, finalvars.c ^ (j/8) * finalvars.alpha / ((finalvars.c ^ (j/8) - 1) * finalvars.alpha + 1)];
    Post_b = [Post_b, finalvars.c ^ (j/8) * (1 - finalvars.alpha) / ((finalvars.c ^ (j/8) - 1) * (1 - finalvars.alpha) + 1)];
end

% drop irrelevant point values
Post_p(to_drop) = [];
Post_b(to_drop) = [];


%% Re-estimate thresholds (really need to do this for the new point values but do it for all to update data structures)
% get thresholds for the regular & bootstrapped posteriors (pathogenic)
finalvars.pthresh{1} = get_all_thresholds(currvars.posteriors_p{i}, thrs, Post_p);

% get thresholds for the regular & bootstrapped posteriors (benign)
finalvars.bthresh{1} = get_all_thresholds(currvars.posteriors_b{i}, flip(thrs), Post_b);

% obtain nonbootstrapped thresholds (pathogenic, benign)
finalvars.ThresholdP{1} = finalvars.pthresh{1}(1, :);
finalvars.ThresholdB{1} = finalvars.bthresh{1}(1, :);

% obtain discounted (one-sided confidence bound-based) thresholds (pathogenic, benign)
finalvars.DiscountedThresholdP{1} = get_discounted_thresholds(finalvars.pthresh{1}(2 : finalvars.B + 1, :), Post_p, finalvars.B, finalvars.discountonesided, 'pathogenic');
finalvars.DiscountedThresholdB{1} = get_discounted_thresholds(finalvars.bthresh{1}(2 : finalvars.B + 1, :), Post_b, finalvars.B, finalvars.discountonesided, 'benign');


%% Update other variables
finalvars.Post_p = Post_p;
finalvars.Post_b = Post_b;
finalvars.j = length(Post_p);

if ~strcmp(orig_file, '')
    finalvars.D = D;
    finalvars.x  = x;
    finalvars.y = y;
    finalvars.g = g;
    finalvars.w = w;
    finalvars.thrs = thrs;
    finalvars.file = file;
    finalvars.files = currvars.files(i); %{orig_file};
    [~, finalvars.filetosave, ~] = fileparts(out_file);
    finalvars.i = i; % only thing that will deviate from the new format so that we can link back to old file-style
    finalvars.increments = currvars.increments(i);
    finalvars.methods = currvars.methods(i);
    finalvars.posteriors_p = currvars.posteriors_p(i);
    finalvars.posteriors_b = currvars.posteriors_b(i);
end


%% Save output
print_thresholds(finalvars.ThresholdP{1}, finalvars.ThresholdB{1}, finalvars.DiscountedThresholdP{1}, finalvars.DiscountedThresholdB{1});
save(out_file, '-struct', 'finalvars', '-v7.3');





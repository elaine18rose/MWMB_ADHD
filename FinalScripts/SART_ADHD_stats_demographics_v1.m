%% Demographic stats
clear all;
close all;
%%
if isempty(findstr(pwd,'thandrillon'))==0
    path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    data_path='/Users/thandrillon/Data/ADHD_MW/EEG/';
    preproc_path='/Users/thandrillon/Data/ADHD_MW/Preproc/';
    path_eeglab='/Users/thandrillon/WorkGit/projects/ext/eeglab/';
    path_ICAlabel='/Users/thandrillon/WorkGit/projects/ext/ICLabel/';
    path_FMINSEARCHBND='/Users/thandrillon/Work/local/FMINSEARCHBND';
    path_ExGauss='/Users/thandrillon/WorkGit/projects/ext/exgauss/';
    save_path = '/Users/thandrillon/WorkGit/projects/inprogress/MWMB_ADHD/tables';

else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/desktop/MATLAB_Functions/fieldtrip/';
    data_path = '/Volumes/Seagate/MWMB_ADHD_SART/EEG/';
    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    path_detectSW = '/Volumes/Seagate/MWMB_ADHD_SART/SW_detection/';
    path_eeglab='/Users/elaine/desktop/MATLAB_Functions/eeglab/';
    path_ICAlabel='/Users/elaine/desktop/MATLAB_Functions/ICLabel/';
    path_ExGauss='/Users/elaine/desktop/MATLAB_Functions/exgauss';
    path_FMINSEARCHBND='/Users/elaine/desktop/MATLAB_Functions/FMINSEARCHBND';
    save_path = '/Users/elaine/desktop/Git_Studies/MWMB_ADHD/tables';
    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(path_LSCPtools))
addpath(path_fieldtrip)
ft_defaults;
addpath(genpath(path_ExGauss))
addpath(genpath(path_FMINSEARCHBND))

demo_table = readtable([save_path filesep 'SART_ADHD_behav_demo_v1.txt']);
demo_table.Group = categorical(demo_table.Group);
demo_table.Sex = categorical(demo_table.Sex);

SW_demo_table = readtable([save_path filesep 'SART_ADHD_SW_demo_v1.txt']);
SW_demo_table.Group = categorical(SW_demo_table.Group);
SW_demo_table.Sex = categorical(SW_demo_table.Sex);
SW_demo_table.SubID = categorical(SW_demo_table.SubID);


probe_demo_table = readtable([save_path filesep 'SART_ADHD_probe_demo_v1.txt']);
probe_demo_table.Group = categorical(probe_demo_table.Group);
probe_demo_table.Sex = categorical(probe_demo_table.Sex);
probe_demo_table.SubID = categorical(probe_demo_table.SubID);
%% Age 
unique_subids = unique(demo_table.SubID);
unique_ages = zeros(length(unique_subids), 1);
% Loop through each unique participant ID to get their age
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    participant_age = demo_table.Age(idx);
    participant_group = demo_table.Group(idx); 
    unique_ages(i) = participant_age(1);  % Store the first (and only) age value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_age = unique_ages(unique_groups == 'A');
Control_age = unique_ages(unique_groups == 'C');


% Test for normality 
% Normal Q-Q plot for ADHD group
normplot(ADHD_age);
title('Normal Q-Q Plot for ADHD');
% Normal Q-Q plot for Control group
normplot(Control_age);
title('Normal Q-Q Plot for Control');

%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_age - mean(ADHD_age)) / std(ADHD_age));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_age - mean(Control_age)) / std(Control_age));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);


% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_age, Control_age);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_age, Control_age);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% Sex
% Loop through each unique participant ID to get their sex
unique_sex = demo_table.Sex((demo_table.Group=='A') | (demo_table.Group=='C')); 
% Get the unique sex for each participant by their SubID 
[~, idx] = unique(demo_table.SubID);  % Find the unique SubID index
% Create a filtered table with unique sex data for each participant
filtered_sex = unique_sex(idx); 

contingency_table = crosstab(demo_table.Group(idx), filtered_sex);
% Display the table
disp('Contingency Table for Sex by Group:');
disp(contingency_table);


% Perform the Chi-square test on the contingency table
[chi2_stat, p_val] = chi2gof(contingency_table(:));  % Flatten the contingency table into a vector

% Display the result
disp(['Chi-squared Statistic: ', num2str(chi2_stat)]);
disp(['p-value: ', num2str(p_val)]);

%% Education level (in years)
% Loop through each unique participant ID to get their Edu level
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    participant_edu = demo_table.Education_Years(idx);
    participant_group = demo_table.Group(idx); 
    unique_edu(i) = participant_edu(1);  % Store the first (and only) edu value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_edu = unique_edu(unique_groups == 'A');
ADHD_edu = ADHD_edu(~isnan(ADHD_edu));
Control_edu = unique_edu(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_edu - mean(ADHD_edu)) / std(ADHD_edu));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_edu - mean(Control_edu)) / std(Control_edu));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_edu, Control_edu);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_edu, Control_edu);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% ESS 
% Loop through each unique participant ID to get their ESS
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    participant_ESS = demo_table.Epworth(idx);
    participant_group = demo_table.Group(idx); 
    unique_ESS(i) = participant_ESS(1);  % Store the first (and only) ESS value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_ESS = unique_ESS(unique_groups == 'A');
% ADHD_ESS = ADHD_edu(~isnan(ADHD_edu));
Control_ESS = unique_ESS(unique_groups == 'C');
Control_ESS = Control_ESS(~isnan(Control_ESS));  % remove NaNs

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_ESS - mean(ADHD_ESS)) / std(ADHD_ESS));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_ESS - mean(Control_ESS)) / std(Control_ESS));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_ESS, Control_ESS);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_ESS, Control_ESS);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% CAARS - Inattentiveness (Subscale E) 
% Loop through each unique participant ID 
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    ppt_CAARSInt = demo_table.CAARS_Inattentiveness(idx);
    participant_group = demo_table.Group(idx); 
    unique_CAARSInt(i) = ppt_CAARSInt(1);  % Store the first (and only) CAARS value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_CAARSInt = unique_CAARSInt(unique_groups == 'A');
% ADHD_ESS = ADHD_edu(~isnan(ADHD_edu));
Control_CAARSInt = unique_CAARSInt(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_CAARSInt - mean(ADHD_CAARSInt)) / std(ADHD_CAARSInt));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_CAARSInt - mean(Control_CAARSInt)) / std(Control_CAARSInt));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_CAARSInt, Control_CAARSInt);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_CAARSInt, Control_CAARSInt);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end


%% CAARS - Hyperactive-Impulsive (Subscale F) 
% Loop through each unique participant ID 
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    ppt_CAARSHyp = demo_table.CAARS_Hyperactivity(idx);
    participant_group = demo_table.Group(idx); 
    unique_CAARSHyp(i) = ppt_CAARSHyp(1);  % Store the first (and only) CAARS value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_CAARSHyp = unique_CAARSHyp(unique_groups == 'A');
% ADHD_ESS = ADHD_edu(~isnan(ADHD_edu));
Control_CAARSHyp = unique_CAARSHyp(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_CAARSHyp - mean(ADHD_CAARSHyp)) / std(ADHD_CAARSHyp));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_CAARSHyp - mean(Control_CAARSHyp)) / std(Control_CAARSHyp));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_CAARSHyp, Control_CAARSHyp);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_CAARSHyp, Control_CAARSHyp);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% CAARS - Total ADHD Symptoms (Subscale G) 
% Loop through each unique participant ID 
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    ppt_CAARSTot = demo_table.CAARS_TotalADHDSymptoms(idx);
    participant_group = demo_table.Group(idx); 
    unique_CAARSTot(i) = ppt_CAARSTot(1);  % Store the first (and only) CAARS value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_CAARSTot = unique_CAARSTot(unique_groups == 'A');
% ADHD_ESS = ADHD_edu(~isnan(ADHD_edu));
Control_CAARSTot = unique_CAARSTot(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_CAARSTot - mean(ADHD_CAARSTot)) / std(ADHD_CAARSTot));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_CAARSTot - mean(Control_CAARSTot)) / std(Control_CAARSTot));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_CAARSTot, Control_CAARSTot);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_CAARSTot, Control_CAARSTot);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% CAARS - ADHD Index (Subscale H) 
% Loop through each unique participant ID 
for i = 1:length(unique_subids)
    idx = strcmp(demo_table.SubID, unique_subids{i});
    ppt_CAARSInd = demo_table.CAARS_ADHDIndex(idx);
    participant_group = demo_table.Group(idx); 
    unique_CAARSInd(i) = ppt_CAARSInd(1);  % Store the first (and only) CAARS value for the participant
    unique_groups(i) = participant_group(1); % Store the group for the participant
end

ADHD_CAARSInd = unique_CAARSInd(unique_groups == 'A');
% ADHD_ESS = ADHD_edu(~isnan(ADHD_edu));
Control_CAARSInd = unique_CAARSInd(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_CAARSInd - mean(ADHD_CAARSInd)) / std(ADHD_CAARSInd));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_CAARSInd - mean(Control_CAARSInd)) / std(Control_CAARSInd));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_CAARSInd, Control_CAARSInd);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal


% t-test
[h, p] = ttest2(ADHD_CAARSInd, Control_CAARSInd);
% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% Average SWD and ESS (correlation)
unique_subids = unique(SW_demo_table.SubID);
nSubs = length(unique_subs);
mean_SWD = zeros(nSubs,1);
ESS_scores = zeros(nSubs,1);

for i = 1:nSubs
    this_idx = sub_idx == i;
    mean_SWD(i) = mean(SW_demo_table.SW_density(this_idx), 'omitnan');
    ESS_scores(i) = SW_demo_table.Epworth(find(this_idx,1)); % get ESS once per subject
end

%Test for normality
%Kolmogorov-Smirnov Test:
[h_SWD, p_SWD] = kstest((mean_SWD - mean(mean_SWD)) / std(mean_SWD));
fprintf('SWD normality: h = %d, p = %.4f\n', h_SWD, p_SWD); % h = 0 means data doesn't sig deviate from normal dist.
[h_ESS, p_ESS] = kstest((ESS_scores - mean(ESS_scores)) / std(ESS_scores));
fprintf('ESS normality: h = %d, p = %.4f\n', h_ESS, p_ESS);


valid_idx = ~isnan(mean_SWD) & ~isnan(ESS_scores);
[r, p] = corr(mean_SWD(valid_idx), ESS_scores(valid_idx), 'Type', 'Pearson'); 
fprintf('Correlation between mean SWD and ESS: r = %.3f, p = %.4f\n', r, p);

% LME for SWD and ESS
lmeSW = fitlme(SW_demo_table, 'SW_density ~ Epworth + (1|SubID)');

%% Average SWAmp and ESS (correlation)
unique_subids = unique(SW_demo_table.SubID);
nSubs = length(unique_subs);
mean_SWAmp= zeros(nSubs,1);
ESS_scores = zeros(nSubs,1);

for i = 1:nSubs
    this_idx = sub_idx == i;
    mean_SWAmp(i) = mean(SW_demo_table.SW_amplitude(this_idx), 'omitnan');
    ESS_scores(i) = SW_demo_table.Epworth(find(this_idx,1)); % get ESS once per subject
end

%Test for normality
%Kolmogorov-Smirnov Test:
[h_SWAmp, p_SWAmp] = kstest((mean_SWAmp - mean(mean_SWAmp)) / std(mean_SWAmp));
fprintf('SWAmp normality: h = %d, p = %.4f\n', h_SWAmp, p_SWAmp); % h = 0 means data doesn't sig deviate from normal dist.
[h_ESS, p_ESS] = kstest((ESS_scores - mean(ESS_scores)) / std(ESS_scores));
fprintf('ESS normality: h = %d, p = %.4f\n', h_ESS, p_ESS);


valid_idx = ~isnan(mean_SWAmp) & ~isnan(ESS_scores);
[r, p] = corr(mean_SWAmp(valid_idx), ESS_scores(valid_idx), 'Type', 'Pearson'); 
fprintf('Correlation between mean SWAmp and ESS: r = %.3f, p = %.4f\n', r, p);

% LME for SWAmp and ESS
lmeSW = fitlme(SW_demo_table, 'SW_amplitude~ Epworth + (1|SubID)');

%% Sleepiness Ratings and ESS

lmeVigxESS = fitlme(probe_demo_table, 'Vigilance~ Epworth + (1|SubID)');

%% ASRS v1.1
unique_subids = unique(probe_demo_table.SubID); % get unique participant IDs
unique_subids(unique_subids == "C015") = [];
unique_ASRS_total = zeros(length(unique_subids),1);
unique_groups = categorical(repmat("", length(unique_subids),1));

for i = 1:length(unique_subids)
    idx = (probe_demo_table.SubID == unique_subids(i));
    asrs_sum = probe_demo_table.ASRS_Inattentiveness(idx) + probe_demo_table.ASRS_Hyperactivity(idx);
    unique_ASRS_total(i) = asrs_sum(1);
    tempGroup = probe_demo_table.Group(idx);
    unique_groups(i) = tempGroup(1);
end

ADHD_ASRS = unique_ASRS_total(unique_groups == 'A');
Control_ASRS = unique_ASRS_total(unique_groups == 'C');

%Test for normality
%Kolmogorov-Smirnov Test:
[h_adhd, p_adhd] = kstest((ADHD_ASRS - mean(ADHD_ASRS)) / std(ADHD_ASRS));
fprintf('ADHD group normality: h = %d, p = %.4f\n', h_adhd, p_adhd); % h = 0 means data doesn't sig deviate from normal dist.
[h_control, p_control] = kstest((Control_ASRS - mean(Control_ASRS)) / std(Control_ASRS));
fprintf('Control group normality: h = %d, p = %.4f\n', h_control, p_control);

% Test for Equal Variance
% F-test for equality of variances
[h_var, p_var] = vartest2(ADHD_ASRS, Control_ASRS);
fprintf('Equality of variances: h = %d, p = %.4f\n', h_var, p_var);% h = 0 means variance is equal; in this case it was unequal


% t-test
if h_var == 1
    fprintf('Variances differ significantly. Using Welch''s t-test.\n');
    [h, p] = ttest2(ADHD_ASRS, Control_ASRS, 'Vartype', 'unequal');
else
    fprintf('Variances equal. Using standard t-test.\n');
    [h, p] = ttest2(ADHD_ASRS, Control_ASRS);
end

% Display results
if h == 0
    fprintf('No significant difference between ADHD and Control groups. p-value = %.4f\n', p);
else
    fprintf('Significant difference between ADHD and Control groups. p-value = %.4f\n', p);
end

%% Sleepiness Ratings and Mind states
stateLabels = {'ON', 'MW', 'MB', 'DK'};
probe_demo_table.StateCat = categorical(probe_demo_table.State, [1 2 3 4], stateLabels);

sleepLabels = {'Ex Alert', 'Alert', 'Sleepy', 'Ex Sleepy'};
probe_demo_table.SleepinessCat = categorical(probe_demo_table.Vigilance, [1 2 3 4], sleepLabels);

[contingency_table, chi2_stat, p_val] = crosstab(probe_demo_table.SleepinessCat, probe_demo_table.StateCat);
disp(['Chi-squared Statistic: ', num2str(chi2_stat)]);
disp(['p-value: ', num2str(p_val)]);
disp(array2table(contingency_table, ...
    'VariableNames', categories(probe_demo_table.StateCat), ...
    'RowNames', categories(probe_demo_table.SleepinessCat)));


%Post hoc 

n = sum(contingency_table(:));                       % Total number of observations
row_totals = sum(contingency_table, 2);              % Row totals (sleepiness levels)
col_totals = sum(contingency_table, 1);              % Column totals (mind states)
expected = row_totals * col_totals / n;              % Expected frequencies

std_residuals = (contingency_table - expected) ./ sqrt(expected); % standardised residuals

resid_table = array2table(std_residuals, ...
    'VariableNames', categories(probe_demo_table.StateCat), ...
    'RowNames', categories(probe_demo_table.SleepinessCat));

disp('Standardized Residuals Table:');
disp(resid_table);

%%% Percentage table for manuscript 
[observed, sleep_levels, mind_states] = crosstab(probe_demo_table.SleepinessCat, probe_demo_table.StateCat);
percent_table = 100 * observed ./ sum(observed, 2);  % each row sums to 100%

disp('Percentage of Mind States by Vigilance Level:');
disp(percent_table);

%%% Report for manuscript
N = sum(observed(:)); % Total N (number of included probes)
df = (size(observed,1) - 1) * (size(observed,2) - 1);% Degrees of freedom = (rows - 1) * (columns - 1)
% Chi-squared statistic manually:
expected = sum(observed,2) * sum(observed,1) / N;
chi2_stat = sum(((observed - expected).^2) ./ expected, 'all');
p_val = 1 - chi2cdf(chi2_stat, df); % p-value

% Calculate Cramér's V
k = min(size(observed)); 
cramers_V = sqrt(chi2_stat / (N * (k - 1)));

% Display everything
fprintf('Chi-squared test of independence:\n');
fprintf('χ²(%d, N = %d) = %.2f, p < %.3g, Cramér''s V = %.3f\n', df, N, chi2_stat, p_val, cramers_V);

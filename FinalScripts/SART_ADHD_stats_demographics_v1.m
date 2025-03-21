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
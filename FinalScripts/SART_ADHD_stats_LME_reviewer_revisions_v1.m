%%  stats
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
    path_fdrbh = '/Users/elaine/desktop/MATLAB_Functions/fdr_bh/';
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
addpath(genpath(path_fdrbh))

demo_table = readtable([save_path filesep 'SART_ADHD_behav_demo_v1.txt']);
demo_table.Group = categorical(demo_table.Group);
demo_table.Sex = categorical(demo_table.Sex);
demo_table.SubID = categorical(demo_table.SubID);
demo_table.Group = reordercats(demo_table.Group, {'C', 'A'});
demo_table(demo_table.SubID == 'C015', :) = []; %participant to exclude

SW_demo_table = readtable([save_path filesep 'SART_ADHD_SW_demo_v1.txt']);
SW_demo_table.Group = categorical(SW_demo_table.Group);
SW_demo_table.Sex = categorical(SW_demo_table.Sex);
SW_demo_table.SubID = categorical(SW_demo_table.SubID);
SW_demo_table(SW_demo_table.SubID == 'C015', :) = []; %participant to exclude


probe_demo_table = readtable([save_path filesep 'SART_ADHD_probe_demo_v1.txt']);
probe_demo_table.Group = categorical(probe_demo_table.Group);
probe_demo_table.Group = reordercats(probe_demo_table.Group, {'C', 'A'});
probe_demo_table.Sex = categorical(probe_demo_table.Sex);
probe_demo_table.SubID = categorical(probe_demo_table.SubID);
probe_demo_table(probe_demo_table.SubID == 'C015', :) = []; %participant to exclude


block_demo_table = readtable([save_path filesep 'SART_ADHD_block_demo_v1.txt']);
block_demo_table.Group = categorical(block_demo_table.Group);
block_demo_table.Sex = categorical(block_demo_table.Sex);
block_demo_table.SubID = categorical(block_demo_table.SubID);
block_demo_table.Group = reordercats(block_demo_table.Group, {'C', 'A'});
block_demo_table(block_demo_table.SubID == 'C015', :) = []; %participant to exclude


%% 20/10/25 Addressing reviewer comments
%FA
mdlFA0 = fitlme(demo_table,'FA~1+BlockN+Group+(1|SubID)');
mdlFA1 = fitlme(demo_table,'FA~1+BlockN*Group+(1|SubID)');
compare(mdlFA0, mdlFA1) % winning model is mdlFA0
mdlFA0_age = fitlme(demo_table,'FA~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlFA0, mdlFA0_age) % winning model is mdlFA0_age
anova(mdlFA0_age)
% mdlFA0_ageInt = fitlme(demo_table,'FA~1+BlockN+Group*Age+(1|SubID)'); % adding age as an interaction with group
% compare(mdlFA0_age, mdlFA0_ageInt)
% anova(mdlFA0_ageInt)

%Misses
mdlmiss0 = fitlme(demo_table,'Misses~1+BlockN+Group+(1|SubID)');
mdlmiss1 = fitlme(demo_table,'Misses~1+BlockN*Group+(1|SubID)');
compare(mdlmiss0, mdlmiss1) % winning model is mdlmiss1
mdlmiss1_age = fitlme(demo_table,'Misses~1+BlockN*Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlmiss1, mdlmiss1_age) % winning model is mdlmiss1_age
anova(mdlmiss1_age)
% mdlmiss1_ageInt = fitlme(demo_table,'Misses~1+BlockN*Group*Age+(1|SubID)'); % adding age as an interaction with group
% compare(mdlmiss1_age, mdlmiss1_ageInt) % winning model is mdlmiss1_ageInt
% anova(mdlmiss1_ageInt)

%RT
mdlRT0  = fitlme(demo_table,'RT~1+BlockN+Group+(1|SubID)'); 
mdlRT1  = fitlme(demo_table,'RT~1+BlockN*Group+(1|SubID)'); 
compare(mdlRT0, mdlRT1) % winning model is mdlRT0
mdlRT0_age = fitlme(demo_table,'RT~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlRT0, mdlRT0_age) % winning model is mdlRT0
anova(mdlRT0)

% RT variability (CV)
mdlcvRT0  = fitlme(block_demo_table,'cvRT~1+BlockN+Group+(1|SubID)'); 
mdlcvRT1  = fitlme(block_demo_table,'cvRT~1+BlockN*Group+(1|SubID)'); 
compare(mdlcvRT0, mdlcvRT1) % winning model is mdlcvRT0
mdlcvRT0_age = fitlme(block_demo_table,'cvRT~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlcvRT0, mdlcvRT0_age) % winning model is mdlcvRT0
anova(mdlcvRT0)


% D prime
mdldprime0 = fitlme(block_demo_table,'dprime~1+BlockN+Group+(1|SubID)');
mdldprime1 = fitlme(block_demo_table,'dprime~1+BlockN*Group+(1|SubID)');  
compare(mdldprime0,mdldprime1) % winning model is mdldprime0
mdldprime0_age = fitlme(block_demo_table,'dprime~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdldprime0,mdldprime0_age) % winning model is mdldprime0
anova(mdldprime0)

%criterion
mdlcrit0 = fitlme(block_demo_table,'criterion~1+BlockN+Group+(1|SubID)'); 
mdlcrit1 = fitlme(block_demo_table,'criterion~1+BlockN*Group+(1|SubID)');  
compare(mdlcrit0,mdlcrit1) % winning model is mdlcrit0
mdlcrit0_age = fitlme(block_demo_table,'criterion~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlcrit0,mdlcrit0_age) % winning model is mdlcrit0
anova(mdlcrit0)

%% Benjamini correction for behaviour p-values 
all_pvals = [ ...
    mdlmiss1_age.Coefficients.pValue(strcmp(mdlmiss1_age.Coefficients.Name, 'Group_A')); ...
    mdlmiss1_age.Coefficients.pValue(strcmp(mdlmiss1_age.Coefficients.Name, 'BlockN')); ...
    mdlmiss1_age.Coefficients.pValue(strcmp(mdlmiss1_age.Coefficients.Name, 'Age')); ...
    mdlmiss1_age.Coefficients.pValue(strcmp(mdlmiss1_age.Coefficients.Name, 'Group_A:BlockN')); ...
    mdlFA0_age.Coefficients.pValue(strcmp(mdlFA0_age.Coefficients.Name, 'Group_A')); ...
    mdlFA0_age.Coefficients.pValue(strcmp(mdlFA0_age.Coefficients.Name, 'BlockN')); ...
    mdlFA0_age.Coefficients.pValue(strcmp(mdlFA0_age.Coefficients.Name, 'Age')); ...
    mdldprime0.Coefficients.pValue(strcmp(mdldprime0.Coefficients.Name, 'Group_A')); ...
    mdldprime0.Coefficients.pValue(strcmp(mdldprime0.Coefficients.Name, 'BlockN')); ...
    mdlcrit0.Coefficients.pValue(strcmp(mdlcrit0.Coefficients.Name, 'Group_A')); ...
    mdlcrit0.Coefficients.pValue(strcmp(mdlcrit0.Coefficients.Name, 'BlockN')); ...
    mdlRT0.Coefficients.pValue(strcmp(mdlRT0.Coefficients.Name, 'Group_A')); ...
    mdlRT0.Coefficients.pValue(strcmp(mdlRT0.Coefficients.Name, 'BlockN')); ...
%     mdlstdRT0.Coefficients.pValue(strcmp(mdlstdRT0.Coefficients.Name, 'Group_A')); ...
%     mdlstdRT0.Coefficients.pValue(strcmp(mdlstdRT0.Coefficients.Name, 'BlockN')); ...
    mdlcvRT0.Coefficients.pValue(strcmp(mdlcvRT0.Coefficients.Name, 'Group_A')); ...
    mdlcvRT0.Coefficients.pValue(strcmp(mdlcvRT0.Coefficients.Name, 'BlockN')) ...
]';

[h, crit_p, ~,adj_p] = fdr_bh(all_pvals, 0.05, 'pdep', 'yes');

labels = {'Miss: Group_A';'Miss: BlockN';'Miss: Age';'Miss: Group_A*BlockN';'FA: Group_A';'FA: BlockN';'FA: Age';'dprime: Group_A';'dprime: BlockN';
    'criterion: Group_A';'criterion: BlockN';'RT: Group_A';'RT: BlockN';
%     'stdRT: Group_A';'stdRT: BlockN';
    'cvRT: Group_A';'cvRT: BlockN'};
behaviour_corr_pV = table(labels, all_pvals(:), adj_p(:), 'VariableNames', ...
    {'Coefficient', 'Raw_p', 'Adj_p'});
behaviour_corr_pV.Significant = behaviour_corr_pV.Adj_p < 0.05;
disp(behaviour_corr_pV)

%%
%%%% mind states %%%%
% On Task
mdlON0 = fitlme(block_demo_table,'ON~1+BlockN+Group+(1|SubID)'); 
mdlON1 = fitlme(block_demo_table,'ON~1+BlockN*Group+(1|SubID)');
compare(mdlON0, mdlON1) % winning model is mdlON1
mdlON0_age = fitlme(block_demo_table,'ON~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlON0, mdlON0_age) % winning model is mdlON0
anova(mdlON0)

% Mind Wandering
mdlMW0 = fitlme(block_demo_table,'MW~1+BlockN+Group+(1|SubID)');
mdlMW1 = fitlme(block_demo_table,'MW~1+BlockN*Group+(1|SubID)');
compare(mdlMW0, mdlMW1)
mdlMW0_age = fitlme(block_demo_table,'MW~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlMW0, mdlMW0_age) % winning model is mdlMW0
anova(mdlMW0)

% Mind Blanking
mdlMB0 = fitlme(block_demo_table,'MB~1+BlockN+Group+(1|SubID)');
mdlMB1 = fitlme(block_demo_table,'MB~1+BlockN*Group+(1|SubID)');
compare(mdlMB0, mdlMB1)
mdlMB0_age = fitlme(block_demo_table,'MB~1+BlockN+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlMB0, mdlMB0_age) % winning model is mdlMB0
anova(mdlMB0)


% Don't Remember
mdlDK0 = fitlme(block_demo_table,'DK~1+BlockN+Group+(1|SubID)');
mdlDK1 = fitlme(block_demo_table,'DK~1+BlockN*Group+(1|SubID)');
compare(mdlDK0, mdlDK1) %NOTE: mdlDK1 was better fitting but the interaction, although sig (p=.02), doesn't withhold sig after bonferroni correction (new alpha = .0125)
mdlDK1_age = fitlme(block_demo_table,'DK~1+BlockN*Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlDK1, mdlDK1_age)
anova(mdlDK1)

%% Benjamini correction for mental states p-values 
all_mental_pvals = [ ...
    mdlON0.Coefficients.pValue(strcmp(mdlON0.Coefficients.Name, 'Group_A')); ...
    mdlON0.Coefficients.pValue(strcmp(mdlON0.Coefficients.Name, 'BlockN')); ...
    mdlMW0.Coefficients.pValue(strcmp(mdlMW0.Coefficients.Name, 'Group_A')); ...
    mdlMW0.Coefficients.pValue(strcmp(mdlMW0.Coefficients.Name, 'BlockN')); ...
    mdlMB0.Coefficients.pValue(strcmp(mdlMB0.Coefficients.Name, 'Group_A')); ...
    mdlMB0.Coefficients.pValue(strcmp(mdlMB0.Coefficients.Name, 'BlockN')); ...
    mdlDK1.Coefficients.pValue(strcmp(mdlDK1.Coefficients.Name, 'Group_A')); ...
    mdlDK1.Coefficients.pValue(strcmp(mdlDK1.Coefficients.Name, 'BlockN')); ...
    mdlDK1.Coefficients.pValue(strcmp(mdlDK1.Coefficients.Name, 'Group_A:BlockN')) ...
]';

[h, crit_p, ~,adj_p] = fdr_bh(all_mental_pvals, 0.05, 'pdep', 'yes');

labels = {'On: Group_A';'On: BlockN';'MW: Group_A';'MW: BlockN';'MB: Group_A';'MB: BlockN';
    'DK: Group_A';'DK: BlockN';'DK: Group_A:BlockN'};
mental_corr_pV = table(labels, all_mental_pvals(:), adj_p(:), 'VariableNames', ...
    {'Coefficient', 'Raw_p', 'Adj_p'});
mental_corr_pV.Significant = mental_corr_pV.Adj_p < 0.05;
disp(mental_corr_pV)

%% 
%%% Sleepiness %%%
mdlvig0 = fitlme(probe_demo_table,'Vigilance~1+Block+Group+(1|SubID)');
mdlvig1 = fitlme(probe_demo_table,'Vigilance~1+Block*Group+(1|SubID)');
compare(mdlvig0, mdlvig1)
mdlvig0_age = fitlme(probe_demo_table,'Vigilance~1+Block+Group+Age+(1|SubID)'); % adding age as a covariate
compare(mdlvig0, mdlvig0_age) % winning model is mdlvig0
% mdlvig0_ageInt = fitlme(probe_demo_table,'Vigilance~1+Block+Group*Age+(1|SubID)'); % adding age as an interaction with group 
% compare(mdlvig0_age, mdlvig0_ageInt) % winning model is mdlvig0_age
anova(mdlvig0)
   



%% 11//11/25 Addressing Reviewer comments about Questionnaire data x SW

results_table = table();

VOI = {'ASRS_Inattentiveness', 'ASRS_Hyperactivity', 'ASRS_5', 'CAARS_Inattentiveness', 'CAARS_Hyperactivity', 'CAARS_TotalADHDSymptoms', 'CAARS_ADHDIndex',... 
    'Epworth', 'MEQ', 'AUDIT', 'DUDIT', 'IPI', 'WURS', 'IDASDysphoria', 'IDASWellBeing', 'IDASPanic', 'IDASDepression'};

for nV = 1:length(VOI)
    % This pools across all electrodes
    mdl = fitlme(SW_demo_table, sprintf('SW_density ~ 1 + %s + Age + Block + (1|SubID) + (1|Elec)', VOI{nV}));
    
    idx_dens = find_trials(mdl.Coefficients.Name, VOI{nV});

    % Extract results
    results_table.Questionnaire{nV} = VOI{nV};
    results_table.Dens_Beta(nV) = mdl.Coefficients.Estimate(idx_dens);
    results_table.Dens_t(nV) = mdl.Coefficients.tStat(idx_dens);
    results_table.Dens_p(nV) = mdl.Coefficients.pValue(idx_dens);
end

[h, crit_p, ~, adj_p] = fdr_bh( results_table.Dens_p, 0.05, 'pdep', 'yes');

results_table.Dens_p_fdr = adj_p;
results_table.Significant = adj_p < 0.05;

results_sorted = sortrows(results_table, 'Dens_p');
disp(results_sorted(:, {'Questionnaire', 'Dens_Beta', 'Dens_t', 'Dens_p', 'Dens_p_fdr', 'Significant'}))


%% p values for Supp Table A
% extracting only 1 datapoint per participant
[~, idx_unique] = unique(SW_demo_table.SubID, 'first');
SW_demo_sub = SW_demo_table(idx_unique, :);

results_q = table();

VOI = {'ASRS_Inattentiveness', 'ASRS_Hyperactivity', 'ASRS_5', 'CAARS_Inattentiveness', 'CAARS_Hyperactivity', ...
       'CAARS_TotalADHDSymptoms', 'CAARS_ADHDIndex', 'Epworth', 'MEQ', 'AUDIT', 'DUDIT', 'IPI', 'WURS', ...
       'IDASDysphoria', 'IDASWellBeing', 'IDASPanic', 'IDASDepression'};

for nV = 1:length(VOI)
    varName = VOI{nV};

    ctrl = SW_demo_sub.(varName)(SW_demo_sub.Group == 'Control');
    adhd = SW_demo_sub.(varName)(SW_demo_sub.Group == 'ADHD');

    % remove NaNs
    ctrl = ctrl(~isnan(ctrl));
    adhd = adhd(~isnan(adhd));

    % Welch’s t-test
    [~, p, ~, stats] = ttest2(ctrl, adhd, 'Vartype', 'unequal');

    % Compute summary stats
    m_ctrl = mean(ctrl); sd_ctrl = std(ctrl);
    m_adhd = mean(adhd); sd_adhd = std(adhd);
    pooledSD = sqrt(((sd_ctrl^2) + (sd_adhd^2)) / 2);
    cohend = (m_adhd - m_ctrl) / pooledSD;

    % Store
    results_q.Questionnaire{nV} = varName;
    results_q.t(nV) = stats.tstat;
    results_q.df(nV) = stats.df;
    results_q.p(nV) = p;
    results_q.mean_CTRL(nV) = m_ctrl;
    results_q.sd_CTRL(nV) = sd_ctrl;
    results_q.mean_ADHD(nV) = m_adhd;
    results_q.sd_CTRL(nV) = sd_ctrl;
    results_q.sd_ADHD(nV) = sd_adhd;
    results_q.sd_ADHD(nV) = sd_adhd;
    results_q.CohensD(nV) = cohend;
end

results_q.Significant = results_q.p < 0.05;

results_q_sorted = sortrows(results_q, 'p');

disp(results_q_sorted(:, {'Questionnaire', 't', 'df', 'p', 'Significant', ...
    'mean_CTRL', 'sd_CTRL', 'mean_ADHD', 'sd_ADHD', 'CohensD'}))


%% percentage of MBs
% Filter only valid rows (in case of NaNs)
valid_rows = ~isnan(probe_demo_table.State);

% Get unique subjects
subs = unique(probe_demo_table.SubID);

% Initialize results table
group_names = {'C', 'A'};
mean_MB = [];
sd_MB = [];

for g = 1:length(group_names)
    group = group_names{g};
    
    % Get participants in this group
    group_subs = unique(probe_demo_table.SubID(strcmp(probe_demo_table.Group, group)));
    
    MB_percent = [];
    for s = 1:length(group_subs)
        sub = group_subs{s};
        sub_data = probe_demo_table(strcmp(probe_demo_table.SubID, sub) & valid_rows, :);
        
        % Calculate %MB for this subject
        total_probes = height(sub_data);
        n_MB = sum(sub_data.State == 3);
        MB_percent(end+1) = (n_MB / total_probes) * 100;
    end
    
    % Store group-level stats
    mean_MB(g) = mean(MB_percent);
    sd_MB(g) = std(MB_percent);
    
    fprintf('%s: Mean MB = %.2f%%, SD = %.2f%%\n', group, mean_MB(g), sd_MB(g));
end




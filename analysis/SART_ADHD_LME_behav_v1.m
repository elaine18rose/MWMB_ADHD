%% Init
clear all; 
close all;

if isempty(findstr(pwd,'thandrillon'))==0 
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    behav_path = '/Users/thandrillon/Data/ADHD_MW/Behaviour/';
    save_path = '/Users/thandrillon/WorkGit/projects/inprogress/MWMB_ADHD/tables';
    path_RainCloudPlot='/Users/thandrillon/WorkGit/projects/ext/RainCloudPlots/';
else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/desktop/MATLAB_Functions/fieldtrip/';
    path_RainCloudPlot='/Users/elaine/desktop/MATLAB_Functions/RainCloudPlots/';
    behav_path = '/Volumes/Seagate/MWMB_ADHD_SART/Behaviour/';
    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    path_detectSW = '/Volumes/Seagate/MWMB_ADHD_SART/SW_detection/';
    save_path = '/Users/elaine/desktop/Git_Studies/MWMB_ADHD/tables';
    
    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
addpath(genpath(path_LSCPtools))
addpath(genpath(path_RainCloudPlot));
addpath(path_fieldtrip)
ft_defaults;
% addpath(genpath(path_ExGauss))
% addpath(genpath(path_FMINSEARCHBND))

addpath(behav_path)

% Loading data 
all_behav_table = readtable([save_path filesep 'MWMB_ADHD_all_behav_byTrial.txt']);
all_probe_table = readtable([save_path filesep 'MWMB_ADHD_all_probe_behav.txt']);
all_block_table = readtable([save_path filesep 'MWMB_ADHD_all_block.txt']);

all_behav_table.Group=categorical(all_behav_table.Group);
all_behav_table.SubID=categorical(all_behav_table.SubID);

all_block_table.Group=categorical(all_block_table.Group);
all_block_table.SubID=categorical(all_block_table.SubID);


%% trial-level performance

% FAs/Commission Errors
mdlFA0  = fitlme(all_behav_table,'FA~1+BlockN+Group+(1|SubID)');
mdlFA1  = fitlme(all_behav_table,'FA~1+BlockN*Group+(1|SubID)'); 
mdlFA2  = fitlme(all_behav_table,'FA~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC and BIC model - Group: p = .02, BlockN p <.001
 %%% Extract fit statistics for each model
AIC_values = [mdlFA0.ModelCriterion.AIC, mdlFA1.ModelCriterion.AIC, mdlFA2.ModelCriterion.AIC];
BIC_values = [mdlFA0.ModelCriterion.BIC, mdlFA1.ModelCriterion.BIC, mdlFA2.ModelCriterion.BIC];
% Display results in a table
ModelNames = {'Model 0', 'Model 1', 'Model 2'};
fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
disp(fit_table);

anova(mdlFA2)


% Misses/Omission Errors
mdlMiss0  = fitlme(all_behav_table,'Misses~1+BlockN+Group+(1|SubID)'); 
mdlMiss1  = fitlme(all_behav_table,'Misses~1+BlockN*Group+(1|SubID)'); 
mdlMiss2  = fitlme(all_behav_table,'Misses~1+BlockN*Group+(BlockN|SubID)'); % Better fitting model - p = .3, BlockN = .0028
% Another way to compare models: 
% compare(mdlMiss0, mdlMiss1)
% compare(mdlMiss1,mdlMiss2)
% compare(mdlMiss0,mdlMiss2)

anova(mdlMiss2)


%RT
mdlRT0  = fitlme(all_behav_table,'RT~1+BlockN+Group+(1|SubID)'); 
mdlRT1  = fitlme(all_behav_table,'RT~1+BlockN*Group+(1|SubID)'); 
mdlRT2  = fitlme(all_behav_table,'RT~1+BlockN*Group+(BlockN|SubID)'); % Winning model; 
 %%% Extract fit statistics for each model
% AIC_values = [mdlRT0.ModelCriterion.AIC, mdlRT1.ModelCriterion.AIC, mdlRT2.ModelCriterion.AIC];
% BIC_values = [mdlRT0.ModelCriterion.BIC, mdlRT1.ModelCriterion.BIC, mdlRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlRT2)


% Standard Deviation of Reaction Times
mdlstdRT0  = fitlme(all_behav_table,'stdRT~1+BlockN+Group+(1|SubID)'); % Winning model
mdlstdRT1  = fitlme(all_behav_table,'stdRT~1+BlockN*Group+(1|SubID)'); 
mdlstdRT2  = fitlme(all_behav_table,'stdRT~1+BlockN*Group+(BlockN|SubID)'); 
% %%% Extract fit statistics for each model
% AIC_values = [mdlstdRT0.ModelCriterion.AIC, mdlstdRT1.ModelCriterion.AIC, mdlstdRT2.ModelCriterion.AIC];
% BIC_values = [mdlstdRT0.ModelCriterion.BIC, mdlstdRT1.ModelCriterion.BIC, mdlstdRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlstdRT0)


% D prime
mdldprime0 = fitlme(all_behav_table,'dprime~1+BlockN+Group+(1|SubID)');
mdldprime1 = fitlme(all_behav_table,'dprime~1+BlockN*Group+(1|SubID)');  % Winning model - Group: p = .037, BlockN: p <.001
mdldprime2 = fitlme(all_behav_table,'dprime~1+BlockN*Group+(BlockN|SubID)');
% %% Extract fit statistics for each model
% AIC_values = [mdldprime0.ModelCriterion.AIC, mdldprime1.ModelCriterion.AIC, mdldprime2.ModelCriterion.AIC];
% BIC_values = [mdldprime0.ModelCriterion.BIC, mdldprime1.ModelCriterion.BIC, mdldprime2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdldprime1)

%criterion
mdlcrit0 = fitlme(all_behav_table,'crit~1+BlockN+Group+(1|SubID)');
mdlcrit1 = fitlme(all_behav_table,'crit~1+BlockN*Group+(1|SubID)');  % Winning model - Group: p = .0, BlockN: p = .999
mdlcrit2 = fitlme(all_behav_table,'crit~1+BlockN*Group+(BlockN|SubID)');
% %% Extract fit statistics for each model
% AIC_values = [mdlcrit0.ModelCriterion.AIC, mdlcrit1.ModelCriterion.AIC, mdlcrit2.ModelCriterion.AIC];
% BIC_values = [mdlcrit0.ModelCriterion.BIC, mdlcrit1.ModelCriterion.BIC, mdlcrit2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlcrit1)

CV_CTR=[];
for nc=1:length(ctrs)
    CV_CTR(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==ctrs(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
CV_ADHD=[];
for nc=1:length(adhds)
    CV_ADHD(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==adhds(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
end

%% block-level stats

% FAs/Commission Errors
mdlblockFA0    = fitlme(all_block_table,'FA~1+BlockN+Group+(1|SubID)'); % Winning BIC model - Group: p = .0086, BlockN: p <.001
mdlblockFA1    = fitlme(all_block_table,'FA~1+BlockN*Group+(1|SubID)');
mdlblockFA2    = fitlme(all_block_table,'FA~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC model - Group: p =.02, BlockN: p <.001
 %%% Extract fit statistics for each model
AIC_values = [mdlblockFA0.ModelCriterion.AIC, mdlblockFA1.ModelCriterion.AIC, mdlblockFA2.ModelCriterion.AIC];
BIC_values = [mdlblockFA0.ModelCriterion.BIC, mdlblockFA1.ModelCriterion.BIC, mdlblockFA2.ModelCriterion.BIC];
% Display results in a table
ModelNames = {'Model 0', 'Model 1', 'Model 2'};
fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
disp(fit_table);

anova(mdlblockFA2)


% Misses/Omission Errors
mdlblockMiss0  = fitlme(all_block_table,'Misses~1+BlockN+Group+(1|SubID)');
mdlblockMiss1  = fitlme(all_block_table,'Misses~1+BlockN*Group+(1|SubID)');
mdlblockMiss2  = fitlme(all_block_table,'Misses~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC and BIC model - Group: p = .34, BlockN: p = .003
 %%% Extract fit statistics for each model
AIC_values = [mdlblockMiss0.ModelCriterion.AIC, mdlblockMiss1.ModelCriterion.AIC, mdlblockMiss2.ModelCriterion.AIC];
BIC_values = [mdlblockMiss0.ModelCriterion.BIC, mdlblockMiss1.ModelCriterion.BIC, mdlblockMiss2.ModelCriterion.BIC];
% Display results in a table
ModelNames = {'Model 0', 'Model 1', 'Model 2'};
fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
disp(fit_table);

anova(mdlblockMiss2)


% Mean RT per block
mdlblockRT0  = fitlme(all_block_table,'RT~1+BlockN+Group+(1|SubID)'); % Winning model for BIC - Group: p = .581, BlockN: p = .324
mdlblockRT1  = fitlme(all_block_table,'RT~1+BlockN*Group+(1|SubID)');
mdlblockRT2  = fitlme(all_block_table,'RT~1+BlockN*Group+(BlockN|SubID)'); % Winning model for AIC - Group: p = .734, BlockN: p = .867
 %%% Extract fit statistics for each model
% AIC_values = [mdlblockRT0.ModelCriterion.AIC, mdlblockRT1.ModelCriterion.AIC, mdlblockRT2.ModelCriterion.AIC];
% BIC_values = [mdlblockRT0.ModelCriterion.BIC, mdlblockRT1.ModelCriterion.BIC, mdlblockRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlblockRT2)


% Standard Deviation of Reaction Times
mdlblockstdRT0 = fitlme(all_block_table,'stdRT~1+BlockN+Group+(1|SubID)'); % Winning model - Group: p = .037, BlockN: p <.001
mdlblockstdRT1 = fitlme(all_block_table,'stdRT~1+BlockN*Group+(1|SubID)');
mdlblockstdRT2 = fitlme(all_block_table,'stdRT~1+BlockN*Group+(BlockN|SubID)');
% %% Extract fit statistics for each model
% AIC_values = [mdlblockstdRT0.ModelCriterion.AIC, mdlblockstdRT1.ModelCriterion.AIC, mdlblockstdRT1.ModelCriterion.AIC];
% BIC_values = [mdlblockstdRT0.ModelCriterion.BIC, mdlblockstdRT1.ModelCriterion.BIC, mdlblockstdRT1.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlblockstdRT0)


% D prime
mdlblockdprime0 = fitlme(all_block_table,'dprime~1+BlockN+Group+(1|SubID)'); % Winning model - Group: p = .037, BlockN: p <.001
mdlblockdprime1 = fitlme(all_block_table,'dprime~1+BlockN*Group+(1|SubID)');
mdlblockdprime2 = fitlme(all_block_table,'dprime~1+BlockN*Group+(BlockN|SubID)');
% %% Extract fit statistics for each model
% AIC_values = [mdlblockstdRT0.ModelCriterion.AIC, mdlblockstdRT1.ModelCriterion.AIC, mdlblockstdRT1.ModelCriterion.AIC];
% BIC_values = [mdlblockstdRT0.ModelCriterion.BIC, mdlblockstdRT1.ModelCriterion.BIC, mdlblockstdRT1.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlblockdprime0)


% Criterion
mdlblockcrit0 = fitlme(all_block_table,'criterion~1+BlockN+Group+(1|SubID)'); % Winning BIC model - Group: p = .0192, BlockN: p = .704
mdlblockcrit1 = fitlme(all_block_table,'criterion~1+BlockN*Group+(1|SubID)');
mdlblockcrit2 = fitlme(all_block_table,'criterion~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC model - Group: p = .0195, BlockN: p = .767
% %% Extract fit statistics for each model
% AIC_values = [mdlblockcrit0.ModelCriterion.AIC, mdlblockcrit1.ModelCriterion.AIC, mdlblockcrit2.ModelCriterion.AIC];
% BIC_values = [mdlblockcrit0.ModelCriterion.BIC, mdlblockcrit1.ModelCriterion.BIC, mdlblockcrit2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlblockcrit2)



% On Task
mdlON0 = fitlme(all_block_table,'ON~1+BlockN+Group+(1|SubID)'); % Winning BIC model; Group: p <.001, BlockN: p <.001
mdlON1 = fitlme(all_block_table,'ON~1+BlockN*Group+(1|SubID)');
mdlON2 = fitlme(all_block_table,'ON~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC model; Group: p <.001, BlockN: p =.005
% %%% Extract fit statistics for each model
% AIC_values = [mdlON0.ModelCriterion.AIC, mdlON1.ModelCriterion.AIC, mdlON2.ModelCriterion.AIC];
% BIC_values = [mdlON0.ModelCriterion.BIC, mdlON1.ModelCriterion.BIC, mdlON2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlON2)


% Mind Wandering
mdlMW0 = fitlme(all_block_table,'MW~1+BlockN+Group+(1|SubID)'); % BIC winning model - Group: p = .0035, BlockN: p <.001
mdlMW1 = fitlme(all_block_table,'MW~1+BlockN*Group+(1|SubID)');
mdlMW2 = fitlme(all_block_table,'MW~1+BlockN*Group+(BlockN|SubID)'); % AIC winning model - Group: p = .0028, BlockN: p = .0074
% %%% Extract fit statistics for each model
% AIC_values = [mdlMW0.ModelCriterion.AIC, mdlMW1.ModelCriterion.AIC, mdlMW2.ModelCriterion.AIC];
% BIC_values = [mdlMW0.ModelCriterion.BIC, mdlMW1.ModelCriterion.BIC, mdlMW2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlMW2)


% Mind Blanking
mdlMB0 = fitlme(all_block_table,'MB~1+BlockN+Group+(1|SubID)');
mdlMB1 = fitlme(all_block_table,'MB~1+BlockN*Group+(1|SubID)');
mdlMB2 = fitlme(all_block_table,'MB~1+BlockN*Group+(BlockN|SubID)'); % Winning model - Group: p = .028, BlockN: p = .1549
%%% Extract fit statistics for each model
% AIC_values = [mdlMB0.ModelCriterion.AIC, mdlMB1.ModelCriterion.AIC, mdlMB2.ModelCriterion.AIC];
% BIC_values = [mdlMB0.ModelCriterion.BIC, mdlMB1.ModelCriterion.BIC, mdlMB2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlMB2)


% Don't Remember
mdlDK0 = fitlme(all_block_table,'DK~1+BlockN+Group+(1|SubID)');
mdlDK1 = fitlme(all_block_table,'DK~1+BlockN*Group+(1|SubID)');
mdlDK2 = fitlme(all_block_table,'DK~1+BlockN*Group+(BlockN|SubID)'); % Winning model; Group: p = 0.83, BlockN: p = .0016
%%% Extract fit statistics for each model
AIC_values = [mdlDK0.ModelCriterion.AIC, mdlDK1.ModelCriterion.AIC, mdlDK2.ModelCriterion.AIC];
BIC_values = [mdlDK0.ModelCriterion.BIC, mdlDK1.ModelCriterion.BIC, mdlDK2.ModelCriterion.BIC];
% Display results in a table
ModelNames = {'Model 0', 'Model 1', 'Model 2'};
fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
disp(fit_table);





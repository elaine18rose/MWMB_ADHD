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
    path_chi2test= '/Users/elaine/desktop/MATLAB_Functions/chi2test/';
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
addpath(path_chi2test)
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

% Loading table with demo + questionnaire + byTrial data 
behav_demo_table = readtable([save_path filesep 'SART_ADHD_behav_demo_v1.txt']);
behav_demo_table.Group=categorical(behav_demo_table.Group);
behav_demo_table.SubID=categorical(behav_demo_table.SubID);
behav_demo_table.Sex=categorical(behav_demo_table.Sex);
behav_demo_table.ReportedSubtype=categorical(behav_demo_table.ReportedSubtype);
behav_demo_table.DIVASubtype=categorical(behav_demo_table.DIVASubtype);

all_probe_table.SubID = string(all_probe_table.SubID); % Convert cell array to string
all_probe_table.Group = string(all_probe_table.Group);

all_probe_table.StateC=categorical(nan(size(all_probe_table,1),1));
all_probe_table.StateC(all_probe_table.State==1)='ON';
all_probe_table.StateC(all_probe_table.State==2)='MW';
all_probe_table.StateC(all_probe_table.State==3)='MB';
all_probe_table.StateC(all_probe_table.State==4)='DK';
all_probe_table.VigC=categorical(nan(size(all_probe_table,1),1));
all_probe_table.VigC(all_probe_table.Vigilance==1)='Ex. Alert';
all_probe_table.VigC(all_probe_table.Vigilance==2)='Alert';
all_probe_table.VigC(all_probe_table.Vigilance==3)='Sleepy';
all_probe_table.VigC(all_probe_table.Vigilance==4)='Ex. Sleepy';
behav_demo_table.SubID = string(behav_demo_table.SubID); % Convert character array to string
behav_demo_summary = varfun(@(x) x(1), behav_demo_table, ...
    'GroupingVariables', 'SubID', ...
    'InputVariables', {'ReportedSubtype'}); %Summary table because there were SubID repeats as it was byTrial
behav_demo_summary.Properties.VariableNames{'Fun_ReportedSubtype'} = 'ReportedSubtype';
probe_subtype_table = join(all_probe_table, behav_demo_summary(:, {'SubID', 'ReportedSubtype'}), 'Keys', 'SubID');
probe_subtype_table.SubID = categorical(probe_subtype_table.SubID); 
probe_subtype_table.Group = categorical(probe_subtype_table.Group); 
%% trial-level performance

% FAs/Commission Errors
mdlFA0  = fitlme(all_behav_table,'FA~1+BlockN+Group+(1|SubID)');
mdlFA1  = fitlme(all_behav_table,'FA~1+BlockN*Group+(1|SubID)'); 
%mdlFA2  = fitlme(all_behav_table,'FA~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC and BIC model - Group: p = .02, BlockN p <.001
 %%% Extract fit statistics for each model
% AIC_values = [mdlFA0.ModelCriterion.AIC, mdlFA1.ModelCriterion.AIC];%, mdlFA2.ModelCriterion.AIC];
% BIC_values = [mdlFA0.ModelCriterion.BIC, mdlFA1.ModelCriterion.BIC];%, mdlFA2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1'};%, 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
compare(mdlFA0, mdlFA1)

anova(mdlFA0)


% Misses/Omission Errors
mdlMiss0  = fitlme(all_behav_table,'Misses~1+BlockN+Group+(1|SubID)'); 
mdlMiss1  = fitlme(all_behav_table,'Misses~1+BlockN*Group+(1|SubID)'); 
% mdlMiss2  = fitlme(all_behav_table,'Misses~1+BlockN*Group+(BlockN|SubID)'); % Better fitting model - p = .3, BlockN = .0028
% Another way to compare models: 
compare(mdlMiss0, mdlMiss1)
% compare(mdlMiss1,mdlMiss2)
% compare(mdlMiss0,mdlMiss2)

anova(mdlMiss1)


%RT
mdlRT0  = fitlme(all_behav_table,'RT~1+BlockN+Group+(1|SubID)'); 
mdlRT1  = fitlme(all_behav_table,'RT~1+BlockN*Group+(1|SubID)'); 
% mdlRT2  = fitlme(all_behav_table,'RT~1+BlockN*Group+(BlockN|SubID)'); % Winning model; 
 %%% Extract fit statistics for each model
% AIC_values = [mdlRT0.ModelCriterion.AIC, mdlRT1.ModelCriterion.AIC, mdlRT2.ModelCriterion.AIC];
% BIC_values = [mdlRT0.ModelCriterion.BIC, mdlRT1.ModelCriterion.BIC, mdlRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
compare(mdlRT0, mdlRT1)

anova(mdlRT0)


% Standard Deviation of Reaction Times; NOTE: Changed to byBlock data 
mdlstdRT0  = fitlme(all_block_table,'stdRT~1+BlockN+Group+(1|SubID)'); 
mdlstdRT1  = fitlme(all_block_table,'stdRT~1+BlockN*Group+(1|SubID)'); 
% mdlstdRT2  = fitlme(all_block_table,'stdRT~1+BlockN*Group+(BlockN|SubID)'); % Winning model - Group = .02, Block = .0058
%%% Extract fit statistics for each model
% AIC_values = [mdlstdRT0.ModelCriterion.AIC, mdlstdRT1.ModelCriterion.AIC, mdlstdRT2.ModelCriterion.AIC];
% BIC_values = [mdlstdRT0.ModelCriterion.BIC, mdlstdRT1.ModelCriterion.BIC, mdlstdRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
compare(mdlstdRT0, mdlstdRT1)

anova(mdlstdRT0)


% D prime
mdldprime0 = fitlme(all_block_table,'dprime~1+BlockN+Group+(1|SubID)');% Winning BIC model - Group: p = .034, BlockN: p <.001
mdldprime1 = fitlme(all_block_table,'dprime~1+BlockN*Group+(1|SubID)');  
mdldprime2 = fitlme(all_block_table,'dprime~1+BlockN*Group+(BlockN|SubID)');% Winning AIC model - Group: p = .256, BlockN: p <.001
% %% Extract fit statistics for each model
% AIC_values = [mdldprime0.ModelCriterion.AIC, mdldprime1.ModelCriterion.AIC, mdldprime2.ModelCriterion.AIC];
% BIC_values = [mdldprime0.ModelCriterion.BIC, mdldprime1.ModelCriterion.BIC, mdldprime2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdldprime0)

%criterion
mdlcrit0 = fitlme(all_block_table,'criterion~1+BlockN+Group+(1|SubID)'); % Winning BIC model - Group: p = .02, BlockN: p = .70
mdlcrit1 = fitlme(all_block_table,'criterion~1+BlockN*Group+(1|SubID)');  
mdlcrit2 = fitlme(all_block_table,'criterion~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC model - Group: p = .02, BlockN: p = .77
% %% Extract fit statistics for each model
% AIC_values = [mdlcrit0.ModelCriterion.AIC, mdlcrit1.ModelCriterion.AIC, mdlcrit2.ModelCriterion.AIC];
% BIC_values = [mdlcrit0.ModelCriterion.BIC, mdlcrit1.ModelCriterion.BIC, mdlcrit2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlcrit1)

% CV_CTR=[];
% for nc=1:length(ctrs)
%     CV_CTR(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==ctrs(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
% end
% CV_ADHD=[];
% for nc=1:length(adhds)
%     CV_ADHD(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==adhds(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
% end

%% block-level stats

% FAs/Commission Errors
% mdlblockFA0    = fitlme(all_block_table,'FA~1+BlockN+Group+(1|SubID)'); % Winning BIC model - Group: p = .0086, BlockN: p <.001
% mdlblockFA1    = fitlme(all_block_table,'FA~1+BlockN*Group+(1|SubID)');
mdlblockFA2    = fitlme(all_block_table,'FA~1+BlockN*Group+(BlockN|SubID)'); % Winning AIC model - Group: p =.02, BlockN: p <.001
 %%% Extract fit statistics for each model
% AIC_values = [mdlblockFA0.ModelCriterion.AIC, mdlblockFA1.ModelCriterion.AIC, mdlblockFA2.ModelCriterion.AIC];
% BIC_values = [mdlblockFA0.ModelCriterion.BIC, mdlblockFA1.ModelCriterion.BIC, mdlblockFA2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

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
compare(mdlON0, mdlON1)

anova(mdlON0)


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
compare(mdlMW0, mdlMW1)

anova(mdlMW0)


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
compare(mdlMB0, mdlMB1)

anova(mdlMB2)


% Don't Remember
mdlDK0 = fitlme(all_block_table,'DK~1+BlockN+Group+(1|SubID)');
mdlDK1 = fitlme(all_block_table,'DK~1+BlockN*Group+(1|SubID)');
mdlDK2 = fitlme(all_block_table,'DK~1+BlockN*Group+(BlockN|SubID)'); % Winning model; Group: p = 0.85, BlockN: p = .0015
% %%% Extract fit statistics for each model
% AIC_values = [mdlDK0.ModelCriterion.AIC, mdlDK1.ModelCriterion.AIC, mdlDK2.ModelCriterion.AIC];
% BIC_values = [mdlDK0.ModelCriterion.BIC, mdlDK1.ModelCriterion.BIC, mdlDK2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
compare(mdlDK0, mdlDK1) %NOTE: mdlDK1 was better fitting but the interaction, although sig (p=.02), doesn't withhold sig after bonferroni correction (new alpha = .0125)

anova(mdlDK1)

%% probe-level stats
all_probe_table.Group=categorical(all_probe_table.Group);
all_probe_table.StateC=categorical(nan(size(all_probe_table,1),1));
all_probe_table.StateC(all_probe_table.State==1)='ON';
all_probe_table.StateC(all_probe_table.State==2)='MW';
all_probe_table.StateC(all_probe_table.State==3)='MB';
all_probe_table.StateC(all_probe_table.State==4)='DK';
all_probe_table.VigC=categorical(nan(size(all_probe_table,1),1));
all_probe_table.VigC(all_probe_table.Vigilance==1)='Ex. Alert';
all_probe_table.VigC(all_probe_table.Vigilance==2)='Alert';
all_probe_table.VigC(all_probe_table.Vigilance==3)='Sleepy';
all_probe_table.VigC(all_probe_table.Vigilance==4)='Ex. Sleepy';

% Intentionality (1 = Entirely intentional to 4 = Entirely unintentional)
% mdlint0 = fitlme(all_probe_table,'Intention~1+StateC+Group+(1|SubID)'); %
% mdlint1 = fitlme(all_block_table,'Intention~1+BlockN*Group+(1|SubID)');
mdlint2 = fitlme(all_probe_table,'Intention~1+StateC*Group+(1|SubID)'); % Winning AIC & BIC model; Group: p =.075; StateC: p <.001; Group*StateC: p <.001
% %%% Extract fit statistics for each model
% AIC_values = [mdlint0.ModelCriterion.AIC, mdlint2.ModelCriterion.AIC];
% %AIC_values = [mdlint0.ModelCriterion.AIC, mdlint1.ModelCriterion.AIC, mdlint2.ModelCriterion.AIC];
% BIC_values = [mdlint0.ModelCriterion.BIC, mdlint2.ModelCriterion.BIC];
% % BIC_values = [mdlint0.ModelCriterion.BIC, mdlint1.ModelCriterion.BIC, mdlint2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1'};
% %ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlint2)


figure;
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

all_probe_table.Group = reordercats(all_probe_table.Group, {'C', 'A'});
% myReportedSubtype=categories(probe_subtype_table.ReportedSubtype);
myGroups=unique(all_probe_table.Group);
myStates=unique(all_probe_table.StateC);

hb=[];
Int_Paired_test=[];
for nSta=1:4
    for_paired_states={};
    for nG=1:2
        % Extract intention values for the current group and state
        state_values = all_probe_table.Intention(all_probe_table.StateC == myStates(nSta) & all_probe_table.Group == myGroups(nG));
        % Calculate percentage for current state and group
        state_percentages(nG, nSta) = numel(state_values) / numel(all_probe_table.Intention(all_probe_table.Group == myGroups(nG))) * 100;
        
        temp_bar=grpstats(all_probe_table.Intention(all_probe_table.StateC==myStates(nSta) & all_probe_table.Group==myGroups(nG)),...
            all_probe_table.SubID(all_probe_table.StateC==myStates(nSta) & all_probe_table.Group==myGroups(nG)));
        hb(nG)=simpleBarPlot(nSta+(2*nG-3)*0.2,temp_bar,Colors(nG,:),0.35,'none',[],2); %none was 'k' before to show SE (or was it SD?) bars
        
        % Add percentage text on top of each bar
        text(nSta + (2 * nG - 3) * 0.2, nanmean(temp_bar) + 0.01, sprintf('%.1f%%', state_percentages(nG, nSta)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'FontWeight', 'bold', 'Color', Colors(nG,:));
        

        for_paired_states{nG}=temp_bar;
    end
    [h,p,ci,stats] = ttest2(for_paired_states{1},for_paired_states{2});
    Int_Paired_test(nSta,1)=p;
    Int_Paired_test(nSta,2)=stats.tstat;
    Int_Paired_test(nSta,3)=stats.df;
end
xlabels = {'On Task', 'Mind Wandering','Mind Blanking', 'Dont Remember'};
ylabels = {'Entirely Intentional', 'Somewhat Intentional','Somewhat Unintentional', 'Entirely Unintentional'};
set(gca,'Xtick',1:4,'XTickLabel',xlabels);
set(gca,'Ytick',1:4,'YTickLabel',ylabels);
ytickangle(45);
xtickangle(45);
legend(hb,myGroups)
format_fig;
ylabel('Intentionality')


% Intention x MW
ctrs=unique(all_behav_table.SubID(all_behav_table.Group=='C' ));
adhds=unique(all_behav_table.SubID(all_behav_table.Group=='A' ));

 numbers_of_interest = [1, 2, 3, 4];
    ADHD_int = [];
    Control_int =[];
    Control_int = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_int = zeros(length(adhds), length(numbers_of_interest));

    clear int_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == ctrs(nc) & probe_subtype_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_int(index, :) = histcounts(int_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_int index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_int_total_occurrences = sum(Control_int);
    % Calculate the percentage distribution
    Ctr_int_percentage_distribution = Ctr_int_total_occurrences / sum(Ctr_int_total_occurrences) * 100;

    clear int_values
    for nc = 1:length(adhds)
        % Extract the State column for the current participant
        int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == adhds(nc) & probe_subtype_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_int(nc, :) = histcounts(int_values, [numbers_of_interest, Inf]);
    end

    ADHD_int_total_occurrences = sum(ADHD_int);
    ADHD_int_percentage_distribution = ADHD_int_total_occurrences / sum(ADHD_int_total_occurrences) * 100;



    % Plot the distribution - bar graphs
    labels = {'Entirely Intentional', 'Somewhat Intentional','Somewhat Unintentional', 'Entirely Unintentional'};
    numGroups = 2;  % Control, ADHD
    groupOffsets = linspace(-0.2, 0.2, numGroups);  % Calculate offsets for each group
    barWidth = 0.29;
    figure%('Position',[100         100         800         600]);

    barPositions = (1:numel(labels)) - 0.3; 
    h1 = bar(barPositions, Ctr_int_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_int_percentage_distribution, 1)
            text(barPositions(i), Ctr_int_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_int_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 25,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    hold on;

    barPositions = (1:numel(labels)) + 0.0;
    h2 = bar(barPositions, ADHD_int_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_int_percentage_distribution, 1)
            text(barPositions(i), ADHD_int_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_int_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(2,:),'FontWeight','bold');
        end
    end

    % Set x-axis labels and ticks
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    xlim([0.5 numel(labels) + 0.5]);
    ylabel('% of Intention Responses for MW');
    ylim([0 80]);
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
    legend([h1, h2], {'Control', 'ADHD'}, 'Location', 'northeast', 'FontSize', 18);
    format_fig;

%% t test for MW x intentionality
ctrs = unique(all_behav_table.SubID(all_behav_table.Group == 'C'));
adhds = unique(all_behav_table.SubID(all_behav_table.Group == 'A'));

numbers_of_interest = [1, 2, 3, 4];
Control_int = zeros(length(ctrs), length(numbers_of_interest));
ADHD_int = zeros(length(adhds), length(numbers_of_interest));

% Populate Control data
index = 1;
for nc = 1:length(ctrs)
    % Extract intentionality values during MW
    int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == ctrs(nc) & probe_subtype_table.StateC == 'MW');
    Control_int(index, :) = histcounts(int_values, [numbers_of_interest, Inf]);
    index = index + 1;
end

% Populate ADHD data
for nc = 1:length(adhds)
    int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == adhds(nc) & probe_subtype_table.StateC == 'MW');
    ADHD_int(nc, :) = histcounts(int_values, [numbers_of_interest, Inf]);
end

% Calculate percentages
Ctr_int_percent = Control_int ./ sum(Control_int, 2) * 100;
ADHD_int_percent = ADHD_int ./ sum(ADHD_int, 2) * 100;


Int_Paired_test = nan(4, 3); % Preallocate for p-value, t-statistic, and df
for nCat = 1:length(numbers_of_interest)
    % Paired t-test for the current category
    [h, p, ci, stats] = ttest2(Ctr_int_percent(:, nCat), ADHD_int_percent(:, nCat));
    Int_Paired_test(nCat, :) = [p, stats.tstat, stats.df];
end

% Convert the results to a table for display
Int_Paired_test_table = array2table(Int_Paired_test, ...
    'VariableNames', {'p', 't-stat', 'df'}, ...
    'RowNames', {'Entirely Int.', 'Somewhat Int.', 'Somewhat Unint.', 'Entirely Unint.'});

disp(Int_Paired_test_table);


% correcting for multiple comparisons - Bonferroni 
alpha = 0.05; % Set your original significance level
num_tests = length(numbers_of_interest); % Number of comparisons
adjusted_alpha = alpha / num_tests; % Adjusted significance level

significant_idx = Int_Paired_test(:, 1) < adjusted_alpha;
disp('Significant results (Bonferroni corrected):');
disp(Int_Paired_test(significant_idx, :));


%% 
% State (1 = On, 2 = MW, 3 = MB)
mdlstate0 = fitlme(all_probe_table,'State~1+Block+Group+(1|SubID)');
mdlstate1 = fitlme(all_probe_table,'State~1+Block*Group+(1|SubID)');
mdlstate2 = fitlme(all_probe_table,'State~1+Block*Group+(Block|SubID)'); % Winning AIC & BIC model; Group: p <.001, Block: p =.032
%%% Extract fit statistics for each model
% AIC_values = [mdlstate0.ModelCriterion.AIC, mdlstate1.ModelCriterion.AIC, mdlstate2.ModelCriterion.AIC];
% BIC_values = [mdlstate0.ModelCriterion.BIC, mdlstate1.ModelCriterion.BIC, mdlstate2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlstate2)


% Distraction (1 = in the room; 2 = personal; 3 = about the task)
% mdldist0 = fitlme(all_probe_table,'Distraction~1+Block+Group+(1|SubID)');
% mdldist1 = fitlme(all_probe_table,'Distraction~1+Block*Group+(1|SubID)');
mdldist2 = fitlme(all_probe_table,'Distraction~1+Block*Group+(Block|SubID)'); % Winning AIC & BIC model; Group: p = .94, Block: p = .99
%%% Extract fit statistics for each model
% AIC_values = [mdldist0.ModelCriterion.AIC, mdldist1.ModelCriterion.AIC, mdldist2.ModelCriterion.AIC];
% BIC_values = [mdldist0.ModelCriterion.BIC, mdldist1.ModelCriterion.BIC, mdldist2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdldist2)

%%%BELOW IS A DRAFT
%%% t-test to determine if there are group differences for reported distraction
ADHD_dis = all_probe_table.Distraction(all_probe_table.Group== 'A');
ADHD_dis = ADHD_dis(~isnan(ADHD_dis));
ctr_dis = all_probe_table.Distraction(all_probe_table.Group=='C');
ctr_dis = ctr_dis(~isnan(ctr_dis));
distraction_types = unique([ADHD_dis; ctr_dis]);

num_types = numel(distraction_types);
ADHD_props = zeros(1, num_types);
ctr_props = zeros(1, num_types);
% Calculate proportions for each distraction type in both groups
for i = 1:num_types
    ADHD_props(i) = sum(ADHD_dis == distraction_types(i)) / numel(ADHD_dis);
    ctr_props(i) = sum(ctr_dis == distraction_types(i)) / numel(ctr_dis);
end

% Perform separate t-tests for each distraction type
p_values = zeros(1, num_types);
for i = 1:num_types
    [~, p] = ttest2(double(ADHD_dis == distraction_types(i)), double(ctr_dis == distraction_types(i)));
    p_values(i) = p;
end

% Correct for multiple comparisons using Bonferroni correction
bonferroni_corrected_p = p_values * num_types; % Multiply by the number of tests
bonferroni_corrected_p(bonferroni_corrected_p > 1) = 1; % p-values cannot exceed 1

% Display results
disp('Distraction Type | Proportion (ADHD) | Proportion (Control) | Raw p-value | Bonferroni Corrected p-value');
for i = 1:num_types
    fprintf('%15d | %16.2f%% | %20.2f%% | %11.4f | %25.4f\n', ...
        distraction_types(i), ADHD_props(i) * 100, ctr_props(i) * 100, p_values(i), bonferroni_corrected_p(i));
end


%% Intentionality (1 = Entirely intentional to 4 = Entirely unintentional)
% mdlint0 = fitlme(all_probe_table,'Intention~1+Block+Group+(1|SubID)');
% mdlint1 = fitlme(all_probe_table,'Intention~1+Block*Group+(1|SubID)');
mdlint2 = fitlme(all_probe_table,'Intention~1+Block*Group+(Block|SubID)'); % Winning AIC & BIC model; Group: p <.001, Block: p <.001
%%% Extract fit statistics for each model
% AIC_values = [mdlint0.ModelCriterion.AIC, mdlint1.ModelCriterion.AIC, mdlint2.ModelCriterion.AIC];
% BIC_values = [mdlint0.ModelCriterion.BIC, mdlint1.ModelCriterion.BIC, mdlint2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlint2)


%% Vigilance (1 = Extremely alert to 4 = Extremely sleepy)
mdlvig0 = fitlme(all_probe_table,'Vigilance~1+Block+Group+(1|SubID)');
mdlvig1 = fitlme(all_probe_table,'Vigilance~1+Block*Group+(1|SubID)');
mdlvig2 = fitlme(all_probe_table,'Vigilance~1+Block*Group+(Block|SubID)'); % Winning AIC & BIC model; Group: p <.001, Block: p <.001
% %%% Extract fit statistics for each model
% AIC_values = [mdlvig0.ModelCriterion.AIC, mdlvig1.ModelCriterion.AIC, mdlvig2.ModelCriterion.AIC];
% BIC_values = [mdlvig0.ModelCriterion.BIC, mdlvig1.ModelCriterion.BIC, mdlvig2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
compare(mdlvig0, mdlvig1)

anova(mdlvig0)


vigilance_levels = 1:4;

% Initialize result arrays
p_values = NaN(length(vigilance_levels), 1);
Q_values = NaN(length(vigilance_levels), 1);

% Loop through each vigilance level
for vigilance_level = vigilance_levels
    % Get counts for ADHD group (for the current vigilance level)
    adhd_vigilance = sum(all_probe_table.Vigilance == vigilance_level & all_probe_table.Group == 'A');
    % Get counts for Control group (for the current vigilance level)
    control_vigilance = sum(all_probe_table.Vigilance == vigilance_level & all_probe_table.Group == 'C');

    % Get counts for ADHD group that did not select the current vigilance level
    adhd_not_vigilance = sum(all_probe_table.Vigilance ~= vigilance_level & all_probe_table.Group == 'A');
    % Get counts for Control group that did not select the current vigilance level
    control_not_vigilance = sum(all_probe_table.Vigilance ~= vigilance_level & all_probe_table.Group == 'C');

    % Create the contingency table (2x2 matrix)
    count_matrix = [adhd_vigilance, adhd_not_vigilance; 
                    control_vigilance, control_not_vigilance];

    % Perform chi-square test using chi2test
    [p, Q] = chi2test(count_matrix);

    % Store results
    p_values(vigilance_level) = p;
    Q_values(vigilance_level) = Q;
    N_values(vigilance_level) = sum(count_matrix, 'all'); % Compute total sample size

    % Display results for this vigilance level
    disp(['Vigilance Level ', num2str(vigilance_level)]);
    disp('Contingency Table:');
    disp(count_matrix);
    disp(['Chi-square test p-value: ', num2str(p)]);
    disp(['Chi-square statistic (Q): ', num2str(Q)]);
    disp(['Total sample size (N): ', num2str(N_values(vigilance_level))]);
    disp('------------------------------');
end

% Display final results for all vigilance levels
disp('Summary of Chi-square test results for all vigilance levels:');
disp(table(vigilance_levels', p_values, Q_values, 'VariableNames', {'VigilanceLevel', 'pValue', 'QStatistic'}));



%% GLMEs for subtype and behaviour (NOTE: no ADHD ppts reported an impulsive/hyperactive subtype)
% Adding 'control' as a subtype for controls 
behav_demo_table1 = behav_demo_table;
behav_demo_table1.Subtype((behav_demo_table1.Group== 'C')) = {'Control'}; 
% Set 'Control' as the reference group
behav_demo_table1.Subtype = reordercats(behav_demo_table1.Subtype, {'Control', 'Combined', 'Inattention'});

% False alarms 
mdlsubtypeFA0 = fitglme(behav_demo_table1, 'FA ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Binomial', 'Link', 'Logit'); % Winning model
% mdlsubtypeFA1 = fitglme(behav_demo_table1, 'FA ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeFA2 = fitglme(behav_demo_table1, 'FA ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');

%%% Extract fit statistics for each model
% AIC_values = [mdlsubtype0.ModelCriterion.AIC, mdlsubtype1.ModelCriterion.AIC, mdlsubtype2.ModelCriterion.AIC];
% BIC_values = [mdlsubtype0.ModelCriterion.BIC, mdlsubtype1.ModelCriterion.BIC, mdlsubtype2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeFA0);
anova(mdlsubtypeFA0)


% Misses
mdlsubtypeMiss0 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Binomial', 'Link', 'Logit'); % Winning model
% mdlsubtypeMiss1 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeMiss2 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');

% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypeMiss0.ModelCriterion.AIC, mdlsubtypeMiss1.ModelCriterion.AIC, mdlsubtypeMiss2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeMiss0.ModelCriterion.BIC, mdlsubtypeMiss1.ModelCriterion.BIC, mdlsubtypeMiss2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeMiss0);
anova(mdlsubtypeMiss0)


% RT
behav_demo_table2 = rmmissing(behav_demo_table1, 'DataVariables', {'RT'}); % Removing data with NAs for RT
% mdlsubtypeRT0 = fitlme(behav_demo_table2, 'RT ~ Subtype + BlockN + (1|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity'); 
mdlsubtypeRT1 = fitlme(behav_demo_table2, 'RT ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity');% Winning model
% mdlsubtypeRT2 = fitlme(behav_demo_table2, 'RT ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity');
% 
% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypeRT0.ModelCriterion.AIC, mdlsubtypeRT1.ModelCriterion.AIC, mdlsubtypeRT2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeRT0.ModelCriterion.BIC, mdlsubtypeRT1.ModelCriterion.BIC, mdlsubtypeRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeRT1);
anova(mdlsubtypeRT1)


% StdRT
mdlsubtypeStdRT0 = fitglme(behav_demo_table1, 'stdRT ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity'); % Winning model
% mdlsubtypeStdRT1 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                  'Distribution', 'Normal', 'Link', 'Identity');
% mdlsubtypeStdRT2 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity');

% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypeStdRT0.ModelCriterion.AIC, mdlsubtypeStdRT1.ModelCriterion.AIC, mdlsubtypeStdRT2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeStdRT0.ModelCriterion.BIC, mdlsubtypeStdRT1.ModelCriterion.BIC, mdlsubtypeStdRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeStdRT0);
anova(mdlsubtypeStdRT0)


% d Prime
mdlsubtypedprime0 = fitglme(behav_demo_table1, 'dprime ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity'); % Winning model
% mdlsubtypedprime1 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                  'Distribution', 'Normal', 'Link', 'Identity');
% mdlsubtypedprime2 = fitglme(behav_demo_table1, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity');
% 
% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypedprime0.ModelCriterion.AIC, mdlsubtypedprime1.ModelCriterion.AIC, mdlsubtypedprime2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypedprime0.ModelCriterion.BIC, mdlsubtypedprime1.ModelCriterion.BIC, mdlsubtypedprime2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypedprime0);
anova(mdlsubtypedprime0)

%% GLMEs just with ADHD
data_ADHD = behav_demo_table(~(behav_demo_table.Group == 'C'),:); 
data_ADHD.SubID=categorical(data_ADHD.SubID);
data_ADHD.Sex=categorical(data_ADHD.Sex);
data_ADHD.Subtype=categorical(data_ADHD.Subtype);

% False alarm
mdlsubtypeFA0 = fitglme(data_ADHD, 'FA ~ Subtype + BlockN + (1|SubID)', ...
               'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeFA1 = fitglme(data_ADHD, 'FA ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeFA2 = fitglme(data_ADHD, 'FA ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');

% % Extract fit statistics for each model
% AIC_values = [mdlsubtypeFA0.ModelCriterion.AIC, mdlsubtypeFA1.ModelCriterion.AIC, mdlsubtypeFA2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeFA0.ModelCriterion.BIC, mdlsubtypeFA1.ModelCriterion.BIC, mdlsubtypeFA2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);


disp(mdlsubtype0);
anova(mdlsubtype0)


% Misses
mdlsubtypeMiss0 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1|SubID)', ...
               'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeMiss1 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');
% mdlsubtypeMiss2 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Binomial', 'Link', 'Logit');

% % Extract fit statistics for each model
% AIC_values = [mdlsubtypeMiss0.ModelCriterion.AIC, mdlsubtypeMiss1.ModelCriterion.AIC, mdlsubtypeMiss2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeMiss0.ModelCriterion.BIC, mdlsubtypeMiss1.ModelCriterion.BIC, mdlsubtypeMiss2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);


disp(mdlsubtypeMiss0);
anova(mdlsubtypeMiss0)


% RT
data_ADHD1 = rmmissing(data_ADHD, 'DataVariables', {'RT'}); % Removing data with NAs for RT
% mdlsubtypeRT0 = fitglme(data_ADHD1, 'RT ~ Subtype + BlockN + (1|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity'); 
mdlsubtypeRT1 = fitglme(data_ADHD1, 'RT ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity');% Winning model (BIC)
% mdlsubtypeRT2 = fitglme(data_ADHD1, 'RT ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity'); % Winning model (AIC)
% 
% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypeRT0.ModelCriterion.AIC, mdlsubtypeRT1.ModelCriterion.AIC, mdlsubtypeRT2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeRT0.ModelCriterion.BIC, mdlsubtypeRT1.ModelCriterion.BIC, mdlsubtypeRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeRT1);
anova(mdlsubtypeRT1)


% StdRT
mdlsubtypeStdRT0 = fitglme(data_ADHD, 'stdRT ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity'); % Winning model
% mdlsubtypeStdRT1 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                  'Distribution', 'Normal', 'Link', 'Identity');
% mdlsubtypeStdRT2 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity');
% 
% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypeStdRT0.ModelCriterion.AIC, mdlsubtypeStdRT1.ModelCriterion.AIC, mdlsubtypeStdRT2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypeStdRT0.ModelCriterion.BIC, mdlsubtypeStdRT1.ModelCriterion.BIC, mdlsubtypeStdRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypeStdRT0);
anova(mdlsubtypeStdRT0)


% d Prime
mdlsubtypedprime0 = fitglme(data_ADHD, 'dprime ~ Subtype + BlockN + (1|SubID)', ...
                'Distribution', 'Normal', 'Link', 'Identity'); % Winning model
% mdlsubtypedprime1 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN|SubID)', ...
%                  'Distribution', 'Normal', 'Link', 'Identity');
% mdlsubtypedprime2 = fitglme(data_ADHD, 'Misses ~ Subtype + BlockN + (1 + BlockN + Subtype|SubID)', ...
%                 'Distribution', 'Normal', 'Link', 'Identity');
% 
% %% Extract fit statistics for each model
% AIC_values = [mdlsubtypedprime0.ModelCriterion.AIC, mdlsubtypedprime1.ModelCriterion.AIC, mdlsubtypedprime2.ModelCriterion.AIC];
% BIC_values = [mdlsubtypedprime0.ModelCriterion.BIC, mdlsubtypedprime1.ModelCriterion.BIC, mdlsubtypedprime2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

disp(mdlsubtypedprime0);
anova(mdlsubtypedprime0)


%% analysis for ADHD subtype and anxiety and depression 
data_ADHD = behav_demo_table(~(behav_demo_table.Group == 'C'),:); 

agg_data_ADHD = groupsummary(data_ADHD, 'SubID', {'max'}, {'Depression', 'Anxiety'}); % Aggregate data because right now we have inflated datapoints for each ppt
[uniqueSubIDs, idx] = unique(data_ADHD.SubID); % Get the first occurrence of each SubID
agg_data_ADHD.ReportedSubtype = data_ADHD.ReportedSubtype(idx);
agg_data_ADHD.DIVASubtype=data_ADHD.DIVASubtype(idx);

agg_data_ADHD.SubID=categorical(agg_data_ADHD.SubID);
agg_data_ADHD.ReportedSubtype=categorical(agg_data_ADHD.ReportedSubtype);
agg_data_ADHD.DIVASubtype=categorical(agg_data_ADHD.DIVASubtype);

subtypes = unique(agg_data_ADHD.ReportedSubtype); % Change to DIVA subtypes once we have everyone's data
comorbidity_stats = zeros(numel(subtypes), 2); % Initialize a matrix for stats

for i = 1:numel(subtypes)
    idx = agg_data_ADHD.ReportedSubtype == subtypes(i); % Logical index for each subtype
    
    % Calculate proportions for depression and anxiety
    comorbidity_stats(i, 1) = mean(agg_data_ADHD.max_Depression(idx)); % Proportion of depression
    comorbidity_stats(i, 2) = mean(agg_data_ADHD.max_Anxiety(idx));    % Proportion of anxiety
end

bar(comorbidity_stats, 'grouped');
xticks(1:numel(subtypes));
xticklabels(subtypes);
ylabel('Proportion');
legend({'Depression', 'Anxiety'});
title('Proportion of Comorbidities by ADHD Subtype');
format_fig

% Depression
[depression_table, subtypeLabels, varLabels] = crosstab(agg_data_ADHD.ReportedSubtype, agg_data_ADHD.max_Depression);

% Perform Fisher's exact test on the contingency table; I have too few observations for Chi test and was getting NaN values
[p, q] = chi2test(depression_table);  
% [p, h] = fishertest(depression_table);
% disp(['Fisher''s exact test p-value: ', num2str(p)]); % No sig difference between depression and subtype (might just also be underpowered)

% Anxiety
[anxiety_table, subtypeLabels, varLabels] = crosstab(agg_data_ADHD.ReportedSubtype, agg_data_ADHD.max_Anxiety);

% Perform Fisher's exact test on the contingency table; I have too few observations for Chi test and was getting NaN values
[p, q] = chi2test(anxiety_table);  
[p, h] = fishertest(anxiety_table);
disp(['Fisher''s exact test p-value: ', num2str(p)]); % Small sig difference between anxiety and subtype 

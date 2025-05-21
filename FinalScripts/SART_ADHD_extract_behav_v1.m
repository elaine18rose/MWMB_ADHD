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
    figures_path = '/Users/elaine/Desktop/Git_Studies/MWMB_ADHD/behaviour/Figures';
    
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

% select relevant files, here baseline blocks
files=dir([behav_path filesep 'wanderIM_behavres_*.mat']);

% NOTE: test_res columns: 
% nblock this_blockcond thiset ntrial this_seq_trial TargetID thisresp stimonset dur_face_no_rand thisresptime  this_nogo this_go

%NOTE: probe_res columns:
% this_probe this_probetime startProbe nblock this_blockcond ntrial probe_responses(:,1)' probe_responses(:,2)' probe_responses(:,3)' probe_responses(:,4)';

%%
behavres_mat=[];
behavgroup_cond=[];
behavSub_ID=[];
stateres_mat=[];
stategroup_cond=[];
resblock_mat=[];
behavres_table=[];

for nF=1:length(files)
    File_Name=files(nF).name;
    if contains(File_Name, "AUTOSAVE") ==1 %Skips autosaves per block
        continue;
    end
    fprintf('... processing %s\n',File_Name);
    SubN=(File_Name(1:end-4));
    load(File_Name)

    if SubN(19)=='C'
        Group='CTR';
    elseif SubN(19)=='A'
        Group='ADHD';
    elseif SubN == 'wanderIM_behavres_001_20Sep2023-1704' | 'wanderIM_behavres_004_03Oct2023-1131' % These participants were missing their ppt group in the name 
        Group = 'CTR';
    else
        Group='missing';
    end
    
    %%% Blocks, MWMB responses and performance
    this_probe = probe_res(:,1);
    probe_block = probe_res(:,4);
    state_resp = probe_res(:, 19); % Note: 1 = On, 2 = MW, 3 = MB
    distraction_resp = probe_res(:,20); % Note: 1 = in the room; 2 = personal; 3 = about the task
    intentional_resp=probe_res(:,21); % Note: 1 = Entirely intentional to 4 = Entirely unintentional
    vigilance_resp=probe_res(:,22); % Note: 1 = Extremely alert to 4 = Extremely sleepy

    BlockN = test_res(:,1);
    nTrial=test_res(:,4);
    go_trials = test_res(:,12); 
    nogo_trials = test_res(:,11); % 1 = reacted (when they shouldn't have); 0 = didn't react (correct resp);
    CR = (test_res(:,11)); 
    FA=(1-CR); % Thomas changed this so now it only counts the number of times participants react to a NoGo
    Hits = (test_res(:,12)); 
    Miss=(1-Hits); %Thomas also changes this so for Hits, it counts responses for Gos and then for Misses it counts how many are missed
    RT_all = test_res(:,10)-test_res(:,8); 
        RT = RT_all; 
        RT(isnan(test_res(:,12)))=NaN;
        RT(RT<0.150)=NaN; warning('removing trials with RT below 150ms')

    % Compiling into tables 
    this_behav=nan(length(test_res),9);
    this_behav(:,3)=BlockN;
    this_behav(:,4)=nTrial;
    this_behav(:,5)=go_trials;
    this_behav(:,6)=nogo_trials;
    this_behav(:,7)=FA;
    this_behav(:,8)=Miss;
    this_behav(:,9)=RT;
    this_behav(:,10)=CR; % Correct Rejection column; i.e. participant didn't respond to 3

    % d prime and criterion
    hit_trials = this_behav((isnan(this_behav(:,5)) == 0), 5); %%% Get the right column from this_behav (column 5) - added isnan
    FA_trials = this_behav((isnan(this_behav(:,7)) == 0), 7);  %%%
    [dprime, crit] = calc_dprime(hit_trials, FA_trials);
    % Display results
    disp(['dprime: ', num2str(dprime)]);
    disp(['Criterion: ', num2str(crit)]);
    % Saving it in correct column
    this_behav(:,11) = dprime;
    this_behav(:,12) = crit;

    % Reaction Time Variability (standard deviation of RTs)
    valid_RT = RT(~isnan(RT));  % Only use non-NaN RT values
    stdRT = std(valid_RT); % Calculate the standard deviation of valid RTs
    % Adding as a new column
    this_behav(:,13) = stdRT;


    this_state2=nan(length(probe_res),11);
    this_state2(:,3)=this_probe;
    this_state2(:,4)=probe_block;
    this_state2(:,5)=state_resp;
    this_state2(:,6)=distraction_resp;
    this_state2(:,7)=intentional_resp;
    this_state2(:,8)=vigilance_resp;
    this_state2(:,9)=nan;
    this_state2(:,10)=nan;
    this_state2(:,11)=nan;

    %     behavres_mat=[behavres_mat; [nF*ones(size(this_behav,1),1) this_behav]];  % This is making a matrix of the variables obtained and put in this_behav
    %     behavgroup_cond=[behavgroup_cond; repmat({Group},size(this_behav,1),1)]; % This is putting the group condition (whether control or ADHD) in a variable for the length of trials

    stateres_mat=[stateres_mat; [nF*ones(size(this_state2,1),1) this_state2]];  % This is making a matrix of the variables obtained and put in this_state
    stategroup_cond=[stategroup_cond; repmat({Group},size(this_state2,1),1)]; % This is putting the group condition (whether control or ADHD) in a variable for the length of trials
    %% Saving the data by trial per participant
    behav_table=array2table(this_behav,...
        'VariableNames',{'SubID','Group','BlockN','nTrial','GoTrials','NoGoTrials','FA','Misses','RT','CR','dprime','crit','stdRT'});
    behav_table.SubID=categorical(behav_table.SubID);
    behav_table.Group=categorical(behav_table.Group);
    if File_Name(19)=='C' || File_Name(19)=='A'
        SubCond=File_Name(19);
        SubID=File_Name(19:22);
        behav_table.SubID=repmat({SubID},size(behav_table,1),1);
        behav_table.Group=repmat({SubCond},size(behav_table,1),1);
        writetable(behav_table,[save_path filesep 'MWMB_ADHD_behav_' File_Name(19:22) '.txt']);
    else %if strcmp(SubN , 'wanderIM_behavres_001_20Sep2023-1704') || strcmp(SubN , 'wanderIM_behavres_004_03Oct2023-1131')
        error('File does not have a group prefix!!')
        %         SubCond='C';
        %         SubID=(['C' File_Name(19:21)]);
        %                behav_table.SubID=repmat({SubID},size(behav_table,1),1);
        %         behav_table.Group=repmat({SubCond},size(behav_table,1),1);
        %         writetable(behav_table,[save_path filesep 'MWMB_ADHD_behav_C' File_Name(19:21) '.txt']);
    end
    behavres_table=[behavres_table ; behav_table];

    probe_table=array2table(this_state2,...
        'VariableNames',{'SubID','Group','ProbeNo','Block','State','Distraction','Intention','Vigilance','Miss','FA','HitRT'});
    probe_table.SubID=categorical(probe_table.SubID);
    probe_table.Group=categorical(probe_table.Group);
    probe_table.SubID=repmat({SubID},size(probe_table,1),1);
    probe_table.Group=repmat({SubCond},size(probe_table,1),1);
    for nP=1:size(probe_table,1)
        ProbeTrialIndex=probe_res(nP,6);
        BlockIndex=probe_res(nP,4);
        go_trials=test_res(test_res(:,1)==BlockIndex & test_res(:,4)<ProbeTrialIndex & test_res(:,5)~=3,:);
        nogo_trials=test_res(test_res(:,1)==BlockIndex & test_res(:,4)<ProbeTrialIndex & test_res(:,5)==3,:);
        go_trials=go_trials(end-17:end,:);
        nogo_trials=nogo_trials(end-1:end,:);
        probe_table.Misses(nP)=nanmean(go_trials(:,end)==0);
        tempRT=go_trials(:,10)-go_trials(:,8); tempRT(tempRT<0.150)=NaN;
        probe_table.RT(nP)=nanmean(tempRT);
        probe_table.FA(nP)=nanmean(nogo_trials(:,end-1)==0);
%         probe_table.stdRT = grpstats(all_RT, test_res(:,1), 'std');% EP added 14/01/25
    end
    if File_Name(19)=='C' || File_Name(19)=='A'
        writetable(probe_table,[save_path filesep 'MWMB_ADHD_probeBehav_' File_Name(19:22) '.txt']);
    elseif strcmp(SubN , 'wanderIM_behavres_001_20Sep2023-1704') || strcmp(SubN , 'wanderIM_behavres_004_03Oct2023-1131')
        writetable(probe_table,[save_path filesep 'MWMB_ADHD_probeBehav_C' File_Name(19:21) '.txt']);
    end

    % by block data
    block_table=zeros(4,18);
    block_table=array2table(block_table,...
        'VariableNames',{'SubID','Group','BlockN','ON','MW','MB','DK','Dist','Perso','Int','Intention','Vigilance','Miss','FA','HitRT','dprime','criterion','stdRT'});
    block_table.SubID=categorical(block_table.SubID);
    block_table.Group=categorical(block_table.Group);

    block_table.BlockN=(1:4)';
    block_table.Misses=1-grpstats(test_res(:,12),test_res(:,1));
    block_table.FA=1-grpstats(test_res(:,11),test_res(:,1));
    all_RT=test_res(:,10)-test_res(:,8);
    all_RT(test_res(:,5)==3)=NaN;
    block_table.RT=grpstats(all_RT,test_res(:,1));

    % D-Prime and Criterion 
    for nbl=1:4
        hit_trials_perB = this_behav((this_behav(:,3) == nbl) & isnan(this_behav(:,5)) == 0, 5); %%% Get the right column from this_behav (column 5) - added isnan
        FA_trials_perB = this_behav((this_behav(:,3) == nbl) & isnan(this_behav(:,7)) == 0, 7);  %%%
        [dprime, crit] = calc_dprime(hit_trials_perB, FA_trials_perB);

        % Display results
        disp(['dprime: ', num2str(dprime)]);
        disp(['Criterion: ', num2str(crit)]);

        % Saving it in correct row
        block_table{block_table{:, 3} == nbl, 16} = dprime;
        block_table{block_table{:, 3} == nbl, 17} = crit;

    end

    % Reaction Time Variability
    block_table.stdRT = grpstats(all_RT, test_res(:,1), 'std');


    block_table.Intention=grpstats(probe_res(:,21),probe_res(:,4));
    block_table.Vigilance=grpstats(probe_res(:,22),probe_res(:,4));

    block_table.ON=grpstats(probe_res(:,19)==1,probe_res(:,4));
    block_table.MW=grpstats(probe_res(:,19)==2,probe_res(:,4));
    block_table.MB=grpstats(probe_res(:,19)==3,probe_res(:,4));
    block_table.DK=grpstats(probe_res(:,19)==4,probe_res(:,4));

    sub_probe_res=probe_res;%(probe_res(:,19)==2,:);
    block_table.Dist=grpstats(sub_probe_res(:,20)==1,sub_probe_res(:,4),'sum')./grpstats(sub_probe_res(:,19)==2,sub_probe_res(:,4),'sum');
    block_table.Perso=grpstats(sub_probe_res(:,20)==2,sub_probe_res(:,4),'sum')./grpstats(sub_probe_res(:,19)==2,sub_probe_res(:,4),'sum');
    block_table.Int=grpstats(sub_probe_res(:,20)==3,sub_probe_res(:,4),'sum')./grpstats(sub_probe_res(:,19)==2,sub_probe_res(:,4),'sum');
    block_table.SubID=repmat({SubID},size(block_table,1),1);
    block_table.Group=repmat({SubCond},size(block_table,1),1);
    if File_Name(19)=='C' || File_Name(19)=='A'
        writetable(block_table,[save_path filesep 'MWMB_ADHD_blockBehav_' File_Name(19:22) '.txt']);
    elseif strcmp(SubN , 'wanderIM_behavres_001_20Sep2023-1704') || strcmp(SubN , 'wanderIM_behavres_004_03Oct2023-1131')
        writetable(block_table,[save_path filesep 'MWMB_ADHD_blockBehav_C' File_Name(19:21) '.txt']);
    end

    if exist('all_behav_table')==0
        all_behav_table=behav_table;
        all_probe_table=probe_table;
        all_block_table=block_table;
    else
        all_behav_table=[all_behav_table ; behav_table];
        all_probe_table=[all_probe_table ; probe_table];
        all_block_table=[all_block_table ; block_table];
    end
end

%% Getting trial numbers by group
% Calculating trials per participant grouped by both SubID and Group
participant_group_trials = grpstats(behavres_table, {'SubID', 'Group'}, 'numel', 'DataVars', 'nTrial');
group_stats = grpstats(participant_group_trials, 'Group', {'mean', 'std'}, 'DataVars', 'numel_nTrial');
% Display the results for each group
fprintf('Average number of trials for Control group (C): %.2f (SD = %.2f)\n', ...
    group_stats.mean_numel_nTrial(strcmp(group_stats.Group, 'C')), ...
    group_stats.std_numel_nTrial(strcmp(group_stats.Group, 'C')));
fprintf('Average number of trials for ADHD group (A): %.2f (SD = %.2f)\n', ...
    group_stats.mean_numel_nTrial(strcmp(group_stats.Group, 'A')), ...
    group_stats.std_numel_nTrial(strcmp(group_stats.Group, 'A')));

disp('Group statistics for trial numbers:');
disp(group_stats);

% fprintf('\nDetailed breakdown by participant:\n');
% disp(sortrows(participant_group_trials, 'Group'));

% t-test to see if there's a sig diff between no. of trials completed
control_trials = participant_group_trials.numel_nTrial(strcmp(participant_group_trials.Group, 'C'));
adhd_trials = participant_group_trials.numel_nTrial(strcmp(participant_group_trials.Group, 'A'));
[h, p, ci, stats] = ttest2(control_trials, adhd_trials);
fprintf('\nT-test results for difference in trial numbers between groups:\n');
fprintf('t(%d) = %.3f, p = %.3f\n', stats.df, stats.tstat, p);

if p < 0.05
    fprintf('There is a significant difference in the number of trials between Control and ADHD groups.\n');
else
    fprintf('There is no significant difference in the number of trials between Control and ADHD groups.\n');
end

%% Saving table of all data 
writetable(all_behav_table,[save_path filesep 'MWMB_ADHD_all_behav_byTrial.txt']);
writetable(all_probe_table, [save_path filesep 'MWMB_ADHD_all_probe_behav.txt']);
writetable(all_block_table, [save_path filesep 'MWMB_ADHD_all_block.txt']);



%% Figures
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;
all_behav_table.Group=categorical(all_behav_table.Group);
all_behav_table.SubID=categorical(all_behav_table.SubID);

ctrs=unique(all_behav_table.SubID(all_behav_table.Group=='C' ));
Miss_CTR=[];
for nc=1:length(ctrs)
    if strcmp(ctrs(nc), 'C015')
        continue; % Skipping C015 'cause ppt didn't understand instructions
    end
    Miss_CTR(nc) = nanmean(all_behav_table.Misses(all_behav_table.SubID == ctrs(nc)));
end
adhds=unique(all_behav_table.SubID(all_behav_table.Group=='A' ));
Miss_ADHD=[];
for nc=1:length(adhds)
    Miss_ADHD(nc)=nanmean(all_behav_table.Misses(all_behav_table.SubID==adhds(nc)));
end

FA_CTR=[];
for nc=1:length(ctrs)
    if strcmp(ctrs(nc), 'C015')
        continue; % Skipping C015 'cause ppt didn't understand instructions
    end
    FA_CTR(nc)=nanmean(all_behav_table.FA(all_behav_table.SubID==ctrs(nc)));
end
FA_ADHD=[];
for nc=1:length(adhds)
    FA_ADHD(nc)=nanmean(all_behav_table.FA(all_behav_table.SubID==adhds(nc)));
end

stdRT_CTR=[];
for nc=1:length(ctrs)
    stdRT_CTR(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==ctrs(nc)));
end
stdRT_ADHD=[];
for nc=1:length(adhds)
    stdRT_ADHD(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==adhds(nc)));
end

CV_CTR=[];
for nc=1:length(ctrs)
    CV_CTR(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==ctrs(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
CV_ADHD=[];
for nc=1:length(adhds)
    CV_ADHD(nc)=nanmean(all_behav_table.stdRT(all_behav_table.SubID==adhds(nc)))./nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
end

Hit_RT_CTR=[];
for nc=1:length(ctrs)
    if strcmp(ctrs(nc), 'C015')
        continue; % Skipping C015 'cause ppt didn't understand instructions
    end
    Hit_RT_CTR(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
Hit_RT_ADHD=[];
for nc=1:length(adhds)
    Hit_RT_ADHD(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
end

dprime_CTR=[];
for nc=1:length(ctrs)
    if strcmp(ctrs(nc), 'C015')
        continue; % Skipping C015 'cause ppt didn't understand instructions
    end
    dprime_CTR(nc)=nanmean(all_behav_table.dprime(all_behav_table.SubID==ctrs(nc))); %NOTE though we don't need to calculate mean here cause it's the same across all trials cause it needs all to calculate
end
dprime_ADHD=[];
for nc=1:length(adhds)
    dprime_ADHD(nc)=nanmean(all_behav_table.dprime(all_behav_table.SubID==adhds(nc)));
end

crit_CTR=[];
for nc=1:length(ctrs)
    if strcmp(ctrs(nc), 'C015')
        continue; % Skipping C015 'cause ppt didn't understand instructions
    end
    crit_CTR(nc)=nanmean(all_behav_table.crit(all_behav_table.SubID==ctrs(nc))); %NOTE though we don't need to calculate mean here cause it's the same across all trials cause it needs all to calculate
end
crit_ADHD=[];
for nc=1:length(adhds)
    crit_ADHD(nc)=nanmean(all_behav_table.crit(all_behav_table.SubID==adhds(nc)));
end

flag_figures = input('Plot figures? (0 for no and 1 for yes) ');
if flag_figures == 1
    %% Mind Wandering by Distraction
    numbers_of_interest = [1, 2, 3];
    ADHD_mw = [];
    Control_mw =[];
    Control_mw = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_mw = zeros(length(adhds), length(numbers_of_interest));

    clear mw_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        mw_values = all_probe_table.Distraction(all_probe_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3) for the current participant
        Control_mw(index, :) = histcounts(mw_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_mw_total_occurrences = sum(Control_mw);
    % Calculate the percentage distribution
    Ctr_mw_percentage_distribution = Ctr_mw_total_occurrences / sum(Ctr_mw_total_occurrences) * 100;

    clear mw_values
    for nc = 1:length(adhds)
        % Extract the State column for the current participant
        mw_values = all_probe_table.Distraction(all_probe_table.SubID == adhds(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_mw(nc, :) = histcounts(mw_values, [numbers_of_interest, Inf]);
    end

    ADHD_mw_total_occurrences = sum(ADHD_mw);
    ADHD_mw_percentage_distribution = ADHD_mw_total_occurrences / sum(ADHD_mw_total_occurrences) * 100;



    % Plot the distribution
    labels = {'Smth in Room', 'Personal','About task'};

    figure('Position',[1955         361         758         413]);
    bar((1:numel(labels))-0.2, Ctr_mw_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',0.38);
    ylabel('% of Types of MW');
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_mw_percentage_distribution, 1)
            text(i-0.2, Ctr_mw_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_mw_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 18,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    format_fig
    hold on;

    bar((1:numel(labels))+0.2, ADHD_mw_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',0.38);
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_mw_percentage_distribution, 1)
            text(i+0.2, ADHD_mw_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_mw_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18,'Color',Colors(2,:),'FontWeight','bold');
        end
    end
    format_fig
    xlim([0.5 3.5])
    ylim([0 80])

    %% Dot plots - Miss/omission errors
     Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

     figure; hold on;
     subjects = [];
     data_to_plot=[];
     meandata_to_plot=[];
     group_labels={'C','A'};
     %all_block_table.Group=categorical(all_block_table.Group);
     for i = 1:4 % number of repetitions
         for j = 1:2 % number of group
             data_to_plot{i, j} = all_block_table.Misses(all_block_table.BlockN == i & strcmp(all_block_table.Group, group_labels{j}));

         end
     end


     for j = 1:2
                  plot((1:4) + (2*j - 3) * 0.1, cellfun(@(x) 100 * nanmean(x), data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
         for i = 1:4
             simpleDotPlot(i + (2*j - 3) * 0.1, 100 * data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end
     set(gca, 'xtick',1:4); %to change y-axis to percentage
     %title(['Omission Errors per Block']);
     ylabel('% of Trials'); xlabel('Block');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     hold on;

     % Create invisible plot objects for the legend (to represent the colored markers)
     h1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', Colors(1,:), 'MarkerEdgeColor', Colors(1,:), 'MarkerSize', 10); % Control
     h2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', Colors(2,:), 'MarkerEdgeColor', Colors(2,:), 'MarkerSize', 10); % ADHD
     % Add the legend with the manually created plot handles
     legend([h1, h2], {'Neurotypical', 'ADHD'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 16, 'Position', [0.25, 0.8, 0, 0]);

     ylim([0 0.04]*100)
     xlim([0.5 4.5])
     saveas(gcf,fullfile(figures_path, 'Fig1_PanelB_MissesPerBlock.svg'))


     for j = 1:2 % number of group
        subjects = unique(all_block_table.SubID(strcmp(all_block_table.Group, group_labels{j}))); % Getting all subjects in this group
         subject_means = [];

         for s = 1:length(subjects)
             subject_mean = mean(all_block_table.Misses(strcmp(all_block_table.Group, group_labels{j}) & ismember(all_block_table.SubID, subjects(s))));
             subject_means = [subject_means; subject_mean];
         end
         data_to_plot_perS{j} = subject_means;
     end
     figure('Position',[2245         400         260         428])
     for j = 1:2 % number of group
         simpleDotPlot((2*j-3)*0.1,100*data_to_plot_perS{j},200,Colors(j,:),1,'k','o',[],3,0,0,0);
     end
hold on;

     ylim([0 0.04]*100)
     set(gca, 'xtick',[-0.1 0.1],'xticklabel',{'NT','ADHD'}); %to change y-axis to percentage
     %title(['Omission']);
     %     ylabel('% of Omission Errors'); xlabel('Block Number');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     saveas(gcf,fullfile(figures_path,'Fig1_PanelB_MissesAvg.svg'))

     %% FA/Comission 
     figure; hold on;
     subjects = [];
     data_to_plot=[];
     meandata_to_plot=[];
     group_labels={'C','A'};
     %all_block_table.Group=categorical(all_block_table.Group);
     for i = 1:4 % number of repetitions
         for j = 1:2 % number of group
             data_to_plot{i, j} = all_block_table.FA(all_block_table.BlockN==i & strcmp(all_block_table.Group, group_labels{j}));
         end
     end


     for j = 1:2
          plot((1:4) + (2*j - 3) * 0.1, cellfun(@(x) 100 * nanmean(x), data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
         for i = 1:4
             simpleDotPlot(i + (2*j - 3) * 0.1, 100 * data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end
     set(gca, 'xtick',1:4); %to change y-axis to percentage
     %title(['Commission Errors per Block']);
     ylabel('% of Trials'); xlabel('Block');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     ylim([0 0.6]*100)
     xlim([0.5 4.5])
     saveas(gcf,fullfile(figures_path,'Fig1_PanelC_FAPerBlock.svg'))


     for j = 1:2 % number of group
         subjects = unique(all_block_table.SubID(strcmp(all_block_table.Group, group_labels{j}))); % Getting all subjects in this group
         subject_means = [];

         for s = 1:length(subjects)
             subject_mean = mean(all_block_table.FA(strcmp(all_block_table.Group, group_labels{j}) & ismember(all_block_table.SubID, subjects(s))));
             subject_means = [subject_means; subject_mean];
         end
         data_to_plot_perS{j} = subject_means;
     end
     figure('Position',[2245         400         260         428])
     for j = 1:2 % number of group
         simpleDotPlot((2*j-3)*0.1,100*data_to_plot_perS{j},200,Colors(j,:),1,'k','o',[],3,0,0,0);
     end
     ylim([0 0.6]*100)
     set(gca, 'xtick',[-0.1 0.1],'xticklabel',{'NT','ADHD'}); %to change y-axis to percentage
     %title(['Commission']);
     %     ylabel('% of Omission Errors'); xlabel('Block Number');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     saveas(gcf,fullfile(figures_path,'Fig1_PanelC_FAAvg.svg'))

     %% RT
     figure; hold on;
     subjects = [];
     data_to_plot=[];
     meandata_to_plot=[];
     group_labels={'C','A'};
     %all_block_table.Group=categorical(all_block_table.Group);
     for i = 1:4 % number of repetitions
         for j = 1:2 % number of group
             data_to_plot{i, j} = all_block_table.RT(all_block_table.BlockN==i & strcmp(all_block_table.Group,group_labels{j}));
         end
     end


     for j = 1:2
          plot((1:4) + (2*j - 3) * 0.1, cellfun(@(x) nanmean(x), data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
         for i = 1:4
             simpleDotPlot(i + (2*j - 3) * 0.1, data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end
     set(gca, 'xtick',1:4); %to change y-axis to percentage
     %title(['Reaction Time across Block']);
     ylabel('Reaction Time (s)'); xlabel('Block');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     ylim([0.35 0.45])
     xlim([0.5 4.5])
     saveas(gcf,fullfile(figures_path,'Fig1_PanelD_RTPerBlock.svg'))


     for j = 1:2 % number of group
         subjects = unique(all_block_table.SubID(strcmp(all_block_table.Group , group_labels{j}))); % Getting all subjects in this group
         subject_means = [];

         for s = 1:length(subjects)
             subject_mean = mean(all_block_table.RT(strcmp(all_block_table.Group , group_labels{j}) & ismember(all_block_table.SubID, subjects(s))));
             subject_means = [subject_means; subject_mean];
         end
         data_to_plot_perS{j} = subject_means;
     end
     figure('Position',[2245         400         260         428])
     for j = 1:2 % number of group
         simpleDotPlot((2*j-3)*0.1,data_to_plot_perS{j},200,Colors(j,:),1,'k','o',[],3,0,0,0);
     end
     ylim([0.35 0.45])
     set(gca, 'xtick',[-0.1 0.1],'xticklabel',{'NT','ADHD'}); %to change y-axis to percentage
     %title(['Reaction Time (s)']);
     %     ylabel('% of Omission Errors'); xlabel('Block Number');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     saveas(gcf,fullfile(figures_path,'Fig1_PanelD_RTAvg.svg'))

          %% stdRT
     figure; hold on;
     subjects = [];
     data_to_plot=[];
     meandata_to_plot=[];
     group_labels={'C','A'};
     %all_block_table.Group=categorical(all_block_table.Group);
     for i = 1:4 % number of repetitions
         for j = 1:2 % number of group
             data_to_plot{i, j} = all_block_table.stdRT(all_block_table.BlockN==i & strcmp(all_block_table.Group,group_labels{j}));
         end
     end


     for j = 1:2
          plot((1:4) + (2*j - 3) * 0.1, cellfun(@(x) nanmean(x), data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
         for i = 1:4
             simpleDotPlot(i + (2*j - 3) * 0.1, data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end
     set(gca, 'xtick',1:4); %to change y-axis to percentage
     %title(['Std Deviation of RT across Block']);
     ylabel('Std Dev RT (s)'); xlabel('Block');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     ylim([0.05 0.2])
     xlim([0.5 4.5])
     saveas(gcf,fullfile(figures_path, 'Fig1_PanelD_stdRTPerBlock.svg'))


     for j = 1:2 % number of group
         subjects = unique(all_block_table.SubID(strcmp(all_block_table.Group , group_labels{j}))); % Getting all subjects in this group
         subject_means = [];

         for s = 1:length(subjects)
             subject_mean = mean(all_block_table.stdRT(strcmp(all_block_table.Group , group_labels{j}) & ismember(all_block_table.SubID, subjects(s))));
             subject_means = [subject_means; subject_mean];
         end
         data_to_plot_perS{j} = subject_means;
     end
     figure('Position',[2245         400         260         428])
     for j = 1:2 % number of group
         simpleDotPlot((2*j-3)*0.1,data_to_plot_perS{j},200,Colors(j,:),1,'k','o',[],3,0,0,0);
     end
     ylim([0.05 0.2])
     set(gca, 'xtick',[-0.1 0.1],'xticklabel',{'NT','ADHD'}); %to change y-axis to percentage
     %title(['Std Deviation', newline,' of RT (s)']);
     %     ylabel('% of Omission Errors'); xlabel('Block Number');
     format_fig;
     set(gca,'FontSize',22,'FontWeight','bold','LineWidth', 1.5);
     saveas(gcf,fullfile(figures_path, 'Fig1_PanelD_stdRTAvg.svg'))

     %% Mind states
         numbers_of_interest = [1, 2, 3, 4];
    ADHD_state = [];
    Control_state =[];
    %Control_state = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_state = zeros(length(adhds), length(numbers_of_interest));

    clear state_values
    index = 1; % Initialize a separate index for Control_state
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        state_values = all_probe_table.State(all_probe_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_state(index, :) = histcounts(state_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_state_total_occurrences = sum(Control_state);
    % Calculate the percentage distribution
    Ctr_state_percentage_distribution = Ctr_state_total_occurrences / sum(Ctr_state_total_occurrences) * 100;


    clear state_values
    for nc = 1:length(adhds)
        % Extract the State column for the current participant
        state_values = all_probe_table.State(all_probe_table.SubID == adhds(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_state(nc, :) = histcounts(state_values, [numbers_of_interest, Inf]);
    end

    ADHD_state_total_occurrences = sum(ADHD_state);
    ADHD_state_percentage_distribution = ADHD_state_total_occurrences / sum(ADHD_state_total_occurrences) * 100;
     numbers_of_interest = [1, 2, 3, 4];
     figure; hold on;
     group_labels = {'Neurotypical', 'ADHD'};
     num_states = length(numbers_of_interest); % 4 states
     data_to_plot = cell(num_states, 2);

     % Prepare data for dot plot
     for i = 1:num_states
         data_to_plot{i, 1} = Control_state(:, i) ./ sum(Control_state, 2) * 100; % Convert to percentage
         data_to_plot{i, 2} = ADHD_state(:, i) ./ sum(ADHD_state, 2) * 100; % Convert to percentage
     end

     % Plot dot plots
     for j = 1:2 % Control (1) and ADHD (2)
         for i = 1:num_states
             simpleDotPlot(i + (2*j - 3) * 0.1, data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
%          % Plot mean lines
%          plot((1:num_states) + (2*j - 3) * 0.1, cellfun(@nanmean, data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
     end

     % Formatting
     labels = {'On','MW','MB ', 'DR'};
     set(gca, 'XTick', 1:num_states, 'XTickLabel', labels);
     xtickangle(45);
     ylabel('% of Mind State');
     title(['Mind State Distribution']);
     format_fig;
     set(gca, 'FontSize', 22, 'FontWeight', 'bold', 'LineWidth', 1.5);
     xlim([0.5 num_states + 0.5]);
     ylim([0 60])
     
     h1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', Colors(1,:), 'MarkerEdgeColor', Colors(1,:), 'MarkerSize', 10); % Control
     h2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', Colors(2,:), 'MarkerEdgeColor', Colors(2,:), 'MarkerSize', 10); % ADHD
     % Add the legend with the manually created plot handles
     legend([h1, h2], {'Neurotypical', 'ADHD'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 16, 'Position', [0.75, 0.72, 0.1, 0.1]);


     % Save figure
     saveas(gcf,fullfile(figures_path, 'Fig2_PanelA_MindStates.svg'));
     %% Intention x MW
     all_probe_table.Group=categorical(all_probe_table.Group);
     all_probe_table.StateC=categorical(nan(size(all_probe_table,1),1));
     all_probe_table.StateC(all_probe_table.State==1)='ON';
     all_probe_table.StateC(all_probe_table.State==2)='MW';
     all_probe_table.StateC(all_probe_table.State==3)='MB';
     all_probe_table.StateC(all_probe_table.State==4)='DK';

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
         int_values = all_probe_table.Intention(all_probe_table.SubID == ctrs(nc) & all_probe_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
         % Count occurrences of each number (1, 2, 3, 4) for the current participant
         Control_int(index, :) = histcounts(int_values, [numbers_of_interest, Inf]);
         index = index + 1; % Increment the Control_int index only when valid data is added
     end

     %     % Calculate the total occurrences of each number across all participants
     %     Ctr_int_total_occurrences = sum(Control_int);
     %     % Calculate the percentage distribution
     %     Ctr_int_percentage_distribution = Ctr_int_total_occurrences / sum(Ctr_int_total_occurrences) * 100;

     clear int_values
     for nc = 1:length(adhds)
         % Extract the State column for the current participant
         int_values = all_probe_table.Intention(all_probe_table.SubID == adhds(nc) & all_probe_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
         % Count occurrences of each number (1, 2, 3, 4) for the current participant
         ADHD_int(nc, :) = histcounts(int_values, [numbers_of_interest, Inf]);
     end

     %     ADHD_int_total_occurrences = sum(ADHD_int);
     %     ADHD_int_percentage_distribution = ADHD_int_total_occurrences / sum(ADHD_int_total_occurrences) * 100;



     % Plot the distribution - dot plots
     numbers_of_interest = [1, 2, 3, 4]; % Intentionality response categories
     figure; hold on;
     group_labels = {'NT', 'ADHD'};
     num_states = length(numbers_of_interest); % 4 intentionality responses

     % Convert raw counts to percentages per participant
     Control_int_percent = Control_int ./ sum(Control_int, 2) * 100;
     ADHD_int_percent = ADHD_int ./ sum(ADHD_int, 2) * 100;

     % Prepare data for dot plot
     data_to_plot = cell(num_states, 2);
     for i = 1:num_states
         data_to_plot{i, 1} = Control_int_percent(:, i);
         data_to_plot{i, 2} = ADHD_int_percent(:, i);
     end

     % Plot dot plots
     for j = 1:2 % Control (1) and ADHD (2)
         % Plot mean lines
         plot((1:num_states) + (2*j - 3) * 0.1, cellfun(@nanmean, data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);
     end
     for j = 1:2 % Control (1) and ADHD (2)
         for i = 1:num_states
             simpleDotPlot(i + (2*j - 3) * 0.1, data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end

     % Formatting
     labels = {'Ent int', 'Some int', 'Some Unint', 'Ent Unint'};
     set(gca, 'XTick', 1:num_states, 'XTickLabel', labels);
     xtickangle(0);
     ylabel(['% of Intention']);
     title(['Intentionality of Mind Wandering']);
     format_fig;
     set(gca, 'FontSize', 22, 'FontWeight', 'bold', 'LineWidth', 1.5);
     xlim([0.5 num_states + 0.5]);
     ylim([0 60]); yticks(0:10:60);

    saveas(gcf,fullfile(figures_path,'Fig2_PanelB_MWxInt.svg'));

     %% Vigilance
     numbers_of_interest = [1, 2, 3, 4];
     ADHD_vig = [];
    Control_vig =[];
    %Control_vig = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_vig = zeros(length(adhds), length(numbers_of_interest));

    clear vig_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        vig_values = all_probe_table.Vigilance(all_probe_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_vig(index, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_vig_total_occurrences = sum(Control_vig);
    % Calculate the percentage distribution
    Ctr_vig_percentage_distribution = Ctr_vig_total_occurrences / sum(Ctr_vig_total_occurrences) * 100;

    clear vig_values
    for nc = 1:length(adhds)
        % Extract the State column for the current participant
        vig_values = all_probe_table.Vigilance(all_probe_table.SubID == adhds(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_vig(nc, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
    end

    ADHD_vig_total_occurrences = sum(ADHD_vig);
    ADHD_vig_percentage_distribution = ADHD_vig_total_occurrences / sum(ADHD_vig_total_occurrences) * 100;

     figure; hold on;
     group_labels = {'NT', 'ADHD'};
     num_states = length(numbers_of_interest); % 4 states
     data_to_plot = cell(num_states, 2);

     % Prepare data for dot plot
     for i = 1:num_states
         data_to_plot{i, 1} = Control_vig(:, i) ./ sum(Control_vig, 2) * 100; % Convert to percentage
         data_to_plot{i, 2} = ADHD_vig(:, i) ./ sum(ADHD_vig, 2) * 100; % Convert to percentage
     end

     % Plot dot plots
     for j = 1:2 % Control (1) and ADHD (2)
         % Plot mean lines
         plot((1:num_states) + (2*j - 3) * 0.1, cellfun(@nanmean, data_to_plot(:, j)), 'Color', Colors(j,:), 'LineWidth', 4);

         for i = 1:num_states
             simpleDotPlot(i + (2*j - 3) * 0.1, data_to_plot{i, j}, 200, Colors(j,:), 1, 'k', 'o', [], 3, 0, 0, 0);
         end
     end

     % Formatting
     labels = {'Ex Alert', 'Alert','Sleepy', 'Ex Sleepy'};

     set(gca, 'XTick', 1:num_states, 'XTickLabel', labels);
     xtickangle(45);
     ylabel('% of Vigilance Ratings');
     title('Vigilance Ratings');
     format_fig;
     set(gca, 'FontSize', 22, 'FontWeight', 'bold', 'LineWidth', 1.5);
     ax = gca; ax.XAxis.FontSize = 18; set(gca, 'FontWeight', 'bold');
     xlim([0.5 num_states + 0.5]);
     yticks(0:10:60);

     % Save figure
     saveas(gcf,fullfile(figures_path, 'Fig2_PanelC_Vig.svg'));

end


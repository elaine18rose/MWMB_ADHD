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

all_behav_table=[];
all_probe_table=[];
all_block_table=[];

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
    
    SubID=File_Name(19:22);
    if strcmp(SubID,'C015')
        continue;
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
    
    % Saving the data by trial per participant
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
        probe_table.Miss(nP)=nanmean(go_trials(:,end)==0);
        tempRT=go_trials(:,10)-go_trials(:,8); tempRT(tempRT<0.150)=NaN;
        probe_table.HitRT(nP)=nanmean(tempRT);
        probe_table.FA(nP)=nanmean(nogo_trials(:,end-1)==0);
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
    block_table.Miss=1-grpstats(test_res(:,12),test_res(:,1));
    block_table.FA=1-grpstats(test_res(:,11),test_res(:,1));
    all_RT=test_res(:,10)-test_res(:,8);
    all_RT(test_res(:,5)==3)=NaN;
    block_table.HitRT=grpstats(all_RT,test_res(:,1));

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
    Miss_CTR(nc) = nanmean(all_behav_table.Misses(all_behav_table.SubID == ctrs(nc)));
end
adhds=unique(all_behav_table.SubID(all_behav_table.Group=='A' ));
Miss_ADHD=[];
for nc=1:length(adhds)
    Miss_ADHD(nc)=nanmean(all_behav_table.Misses(all_behav_table.SubID==adhds(nc)));
end

FA_CTR=[];
for nc=1:length(ctrs)
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
    Hit_RT_CTR(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
Hit_RT_ADHD=[];
for nc=1:length(adhds)
    Hit_RT_ADHD(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
end

dprime_CTR=[];
for nc=1:length(ctrs)
    dprime_CTR(nc)=nanmean(all_behav_table.dprime(all_behav_table.SubID==ctrs(nc))); %NOTE though we don't need to calculat mean here cause it's the same across all trials cause it needs all to calculate
end
dprime_ADHD=[];
for nc=1:length(adhds)
    dprime_ADHD(nc)=nanmean(all_behav_table.dprime(all_behav_table.SubID==adhds(nc)));
end

%%
all_Miss=100*[Miss_CTR , Miss_ADHD]';
all_Group=[zeros(1,length(Miss_CTR)) , ones(1,length(Miss_ADHD))]';

figure('Position',[347   167   441   447]);
h = daviolinplot(all_Miss,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'CTRL','ADHD'},...
    'legend',{'CTRL','ADHD'},'boxwidth',1.5,'scattersize',70);
ylabel('% Misses');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([-0.5 15])

%%
all_FA=100*[FA_CTR , FA_ADHD]';
all_Group=[zeros(1,length(FA_CTR)) , ones(1,length(FA_ADHD))]';

figure('Position',[347   167   441   447]);
h = daviolinplot(all_FA,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'CTRL','ADHD'},...
    'boxwidth',1.5,'scattersize',70);
ylabel('% False Alarms');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([0 100])


%% Mind States
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



    % Plot the distribution
    labels = {'ON', 'MW','MB', '??'};


    all_States=[Control_state ; ADHD_state]./40*100;
    all_Group=[zeros(size(Control_state,1),1) ; ones(size(ADHD_state,1),1)]';

    figure('Position',[347         335        1028         279]);
    h = daviolinplot(all_States,'groups',all_Group,'outsymbol','k+','colors',Colors,...
        'boxcolors','same','scatter',1,'jitter',1,'xtlabels', labels,...
        'boxwidth',1,'scattersize',30);
    ylabel('% Responses');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
    set(gca,'FontSize',10);
    format_fig;
    ylim([0 100])



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


    all_Vigs=[Control_vig ; ADHD_vig]./40*100;

    % Plot the distribution
    labels = {'Alert ++', 'Alert +','Sleepy +', 'Sleepy ++'};

    figure('Position',[347         335        1028         279]);
    h = daviolinplot(all_Vigs,'groups',all_Group,'outsymbol','k+','colors',Colors,...
        'boxcolors','same','scatter',1,'jitter',1,'xtlabels', labels,...
        'boxwidth',1,'scattersize',30);
    ylabel('% Responses');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
    set(gca,'FontSize',10);
    format_fig;
    ylim([0 100])
%%
%%%%% EXAMPLE OF STATS
% block-level stats
all_block_table.Group=categorical(all_block_table.Group);
all_block_table.SubID=categorical(all_block_table.SubID);

mdlblockFA=fitlme(all_block_table,'FA~1+BlockN+Group+(1|SubID)');
mdlblockMiss=fitlme(all_block_table,'Miss~1+BlockN+Group+(1|SubID)');

mdlON=fitlme(all_block_table,'ON~1+BlockN+Group+(1|SubID)');
mdlMW=fitlme(all_block_table,'MW~1+BlockN+Group+(1|SubID)');
mdlMB=fitlme(all_block_table,'MB~1+BlockN+Group+(1|SubID)');
mdlDK=fitlme(all_block_table,'DK~1+BlockN+Group+(1|SubID)');

mdlVig=fitlme(all_block_table,'Vigilance~1+BlockN+Group+(1|SubID)');

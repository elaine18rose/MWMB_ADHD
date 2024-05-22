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
stateres_mat=[];
stategroup_cond=[];
resblock_mat=[];


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
    blockN = test_res(:,1);

    this_probe = probe_res(:,1);
    probe_block = probe_res(:,4);
    state_resp = probe_res(:, 19); % Note: 1 = On, 2 = MW, 3 = MB
    distraction_resp = probe_res(:,20); % Note: 1 = in the room; 2 = personal; 3 = about the task
    intentional_resp=probe_res(:,21); % Note: 1 = Entirely intentional to 4 = Entirely unintentional
    vigilance_resp=probe_res(:,22); % Note: 1 = Extremely alert to 4 = Extremely sleepy

    go_trials = test_res(:,12); %FLAG: Thomas, please check that I'm coding/interpreting this correctly
    nogo_trials = test_res(:,11); % 1 = reacted (when they shouldn't have); 0 = didn't react (correct resp);
    nTrial=test_res(:,4);
    CR = (test_res(:,11)); FA=(1-CR); % Thomas changed this so now it only counts the number of times participants react to a NoGo
    Hits = (test_res(:,12)); Miss=(1-Hits); %Thomas also changes this so for Hits, it counts responses for Gos and then for Misses it counts how many are missed
    RT_all = test_res(:,10)-test_res(:,8); %FLAG: I just copied this from the general code to check the behav data. This doesn't exclude RTs for FAs
    RT = RT_all; RT(isnan(test_res(:,12)))=NaN;
    RT(RT<0.150)=NaN; warning('removing trials with RT below 150ms')
    %correct_RT = RT(Hits); % FLAG to fix

    % Compiling into tables 
    this_behav=nan(length(test_res),9);
    this_behav(:,3)=blockN;
    this_behav(:,4)=nTrial;
    this_behav(:,5)=go_trials;
    this_behav(:,6)=nogo_trials;
    this_behav(:,7)=FA;
    this_behav(:,8)=Miss;
    this_behav(:,9)=RT;

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
% 
%     stateres_mat=[stateres_mat; [nF*ones(size(this_state,1),1) this_state]];  % This is making a matrix of the variables obtained and put in this_state
%     stategroup_cond=[stategroup_cond; repmat({Group},size(this_state,1),1)]; % This is putting the group condition (whether control or ADHD) in a variable for the length of trials
    %% Saving the data by trial per participant
    behav_table=array2table(this_behav,...
        'VariableNames',{'SubID','Group','BlockN','nTrial','GoTrials','NoGoTrials','FA','Misses','RT'});
    behav_table.SubID=categorical(behav_table.SubID);
behav_table.Group=categorical(behav_table.Group);
    if File_Name(19)=='C' || File_Name(19)=='A'
        SubCond=File_Name(19);
        SubID=File_Name(19:22);
        behav_table.SubID=repmat({SubID},size(behav_table,1),1);
        behav_table.Group=repmat({SubCond},size(behav_table,1),1);
        writetable(behav_table,[save_path filesep 'MWMB_ADHD_behav_' File_Name(19:22) '.txt']);
    elseif strcmp(SubN , 'wanderIM_behavres_001_20Sep2023-1704') || strcmp(SubN , 'wanderIM_behavres_004_03Oct2023-1131')
        SubCond='C';
        SubID=(['C' File_Name(19:21)]);
               behav_table.SubID=repmat({SubID},size(behav_table,1),1);
        behav_table.Group=repmat({SubCond},size(behav_table,1),1);
        writetable(behav_table,[save_path filesep 'MWMB_ADHD_behav_C' File_Name(19:21) '.txt']);
    end

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
    block_table=zeros(4,15);
    block_table=array2table(block_table,...
        'VariableNames',{'SubID','Group','Block','ON','MW','MB','DK','Dist','Perso','Int','Intention','Vigilance','Miss','FA','HitRT'});
    block_table.SubID=categorical(block_table.SubID);
    block_table.Group=categorical(block_table.Group);

    block_table.Block=(1:4)';
    block_table.Miss=1-grpstats(test_res(:,12),test_res(:,1));
    block_table.FA=1-grpstats(test_res(:,11),test_res(:,1));
    all_RT=test_res(:,10)-test_res(:,8);
    all_RT(test_res(:,5)==3)=NaN;
    block_table.HitRT=grpstats(all_RT,test_res(:,1));
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

%     %% Making a table of the by trial behavioural results per participant
%     table1=array2table(behavres_mat,'VariableNames',{'SubID','BlockN','TrialN','GoTrials','NoGoTrials','FA','Misses','RT'});
%     table1.Group=behavgroup_cond;
%     table1.SubID=categorical(table1.SubID);
%     table1.Group=categorical(table1.Group);
%     % table1.Group=reordercats(table1.Group,[2,1]);
% 
%     writetable(table1,[save_path filesep 'MWMB_ADHD_behav_bytrial.txt']); %%% FLAG: Need to fix "SubID" because right now it's wrong
% 
% 
%     table2=array2table(stateres_mat,'VariableNames',{'SubID','ProbeN','BlockN','State','Distraction','Intention','Vigilance'});
%     table2.Group=stategroup_cond;
%     table2.SubID=categorical(table2.SubID);
%     table2.Group=categorical(table2.Group);
%     % table1.Group=reordercats(table1.Group,[2,1]);
% 
%     writetable(table2,[save_path filesep 'MWMB_ADHD_state_bytrial.txt']); %%% FLAG: Need to fix "SubID" because right now it's wrong


%%
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;
all_behav_table.Group=categorical(all_behav_table.Group);
all_behav_table.SubID=categorical(all_behav_table.SubID);

ctrs=unique(all_behav_table.SubID(all_behav_table.Group=='C' ));
Miss_CTR=[];
for nc=1:length(ctrs)
    Miss_CTR(nc)=nanmean(all_behav_table.Misses(all_behav_table.SubID==ctrs(nc)));
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

% stdRT_CTR=[];
% for nc=1:length(ctrs)
%     stdRT_CTR(nc)=nanmean(table1.stdRT(table1.SubID==ctrs(nc)))./nanmean(table1.Hit_RT(table2.SubID==ctrs(nc)));
% end
% stdRT_ADHD=[];
% for nc=1:length(adhds)
%     stdRT_ADHD(nc)=nanmean(table1.stdRT(table1.SubID==adhds(nc)))./nanmean(table1.Hit_RT(table2.SubID==adhds(nc)));
% end

Hit_RT_CTR=[];
for nc=1:length(ctrs)
    Hit_RT_CTR(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
Hit_RT_ADHD=[];
for nc=1:length(adhds)
    Hit_RT_ADHD(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==adhds(nc)));
end


%Miss
figure;
h1 = raincloud_plot(100*Miss_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',5,'bound_data',[0 100]);
h2 = raincloud_plot(100*Miss_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',5,'bound_data',[0 100]);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 20], 'YLim',ylim.*[1 1.5]); % set(gca,'XLim', [-10 110], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('MISS'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'}); 
set(gca,'YColor','none') % removes Y-axis 

[h,p,ci,stats] = ttest2(Miss_CTR,Miss_ADHD);
Misseffect = meanEffectSize(Miss_CTR,Miss_ADHD);
figure; gardnerAltmanPlot(Miss_ADHD,Miss_CTR);
%% FA
figure;
h1 = raincloud_plot(100*FA_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',5,'bound_data',[0 100]);
h2 = raincloud_plot(100*FA_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',5,'bound_data',[0 100]);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 100], 'YLim', ylim.*[1 1.5]); %set(gca,'XLim', [0 6], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('FA'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});
set(gca,'YColor','none') % removes Y-axis 

[h,p,ci,stats] = ttest2(FA_CTR,FA_ADHD);
FAeffect = meanEffectSize(FA_CTR,FA_ADHD);
figure; gardnerAltmanPlot(FA_ADHD,FA_CTR);

%%
% %stdRT
% figure;
% h1 = raincloud_plot(stdRT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.02);
% h2 = raincloud_plot(stdRT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.02);
% 
% set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
% set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
% set(gca,'XLim', [0.0 0.4], 'YLim', ylim.*[1 1.5]); %set(gca,'XLim', [0.0 0.3], 'YLim', ylim.*[0.5 1.7]);
% format_fig; title('stdRT/meanRT'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%Hit_RT
figure;
h1 = raincloud_plot(Hit_RT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.04);
h2 = raincloud_plot(Hit_RT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.04);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0.1 0.8], 'YLim', ylim.*[1 1.5]); %set(gca,'XLim', [1.1 1.9], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('RT (s)'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});
set(gca,'YColor','none') % removes Y-axis 


%RTeffect = meanEffectSize(Hit_RT_CTR,Hit_RT_ADHD);
%figure; gardnerAltmanPlot(Hit_RT_ADHD,Hit_RT_CTR)
%%
%% Repeated Measures plot
% %Miss
% data_to_plot=[];
% group_labels={'C','A'};
% for i = 1:4 % number of repetitions
%     for j = 1:2 % number of group
%         data_to_plot{i, j} = all_behav_table.Misses(all_behav_table.BlockN==i & all_behav_table.Group==group_labels{j});
%     end
% end
% 
% figure; hold on;
% h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'XLim', [0.01 0.08]);
% xtix = get(gca, 'xtick'); %to change y-axis to percentage
% set(gca, 'xtick',xtix, 'xtickLabel',xtix*100); %to change y-axis to percentage
% title(['Misses per Block']);
% xlabel('Percentage of Misses'); ylabel('Block Number');
% format_fig;
% 
% %False Alarms
% data_to_plot=[];
% group_labels={'C','A'};
% for i = 1:4 % number of repetitions
%     for j = 1:2 % number of group
%         data_to_plot{i, j} = all_behav_table.FA(all_behav_table.BlockN==i & all_behav_table.Group==group_labels{j});
%     end
% end
% 
% figure; hold on;
% h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'XLim', [-0.015 1]);
% xtix = get(gca, 'xtick'); %to change y-axis to percentage
% set(gca, 'xtick',xtix, 'xtickLabel',xtix*100); %to change y-axis to percentage
% title(['False alarms per Block']);
% xlabel('Percentage of False Alarms'); ylabel('Block Number');
% format_fig;
% 
% %Hit RT
% data_to_plot=[];
% group_labels={'C','A'};
% for i = 1:4 % number of repetitions
%     for j = 1:2 % number of group
%         data_to_plot{i, j} = all_behav_table.RT(all_behav_table.BlockN==i  & all_behav_table.Group==group_labels{j});
%     end
% end
% 
% figure; hold on;
% h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'XLim', [0 .75]);
% title(['Hit Reaction Times per Block']);
% xlabel('Reaction Times (seconds)'); ylabel('Block Number');
% format_fig;

% %stdRT
% data_to_plot=[];
% group_labels={'CTR','ADHD'};
% for i = 1:4 % number of repetitions
%     for j = 1:2 % number of group
%         data_to_plot{i, j} = all_behav_table.stdRT(all_behav_table.BlockN==i & all_behav_table.Group==group_labels{j});
%     end
% end
% 
% figure; hold on;
% h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'XLim', [-0.3 1.6]);
% title(['stdRTs per block']);
% format_fig;

%% plot the inter-probe interval
% plot(diff(probe_res(:,3)))
% % compute the prooportion of mind states
% nanmean(probe_res(:,19)==1) % on
% nanmean(probe_res(:,19)==2) % MW
% nanmean(probe_res(:,19)==3) % MB
% % compute the vigilance ratings
% nanmean(probe_res(:,22)==1)
% nanmean(probe_res(:,22)==2)
% nanmean(probe_res(:,22)==3)
% nanmean(probe_res(:,22)==4)
% % compute the perf on go trials
% nanmean(test_res(:,12))
% % compute the perf on nogo trials
% nanmean(test_res(:,11))
% % plot the distribution of RTs
% all_RT=test_res(:,10)-test_res(:,8);
% figure; hist(all_RT,100)
% 

%% Box plots of mean distribution 
% Misses
% lineWidth = 3;
% figure; 
% subplot(1,2,1)
% boxplot(100*Miss_CTR,'Colors',Colors(1,:), 'symbol', 'o k')
% ylim([0 15])
% ylabel('Misses')
% xlabel ('Control')
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     set(h(j),'LineWidth',lineWidth)
% 
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.5);
% end
% format_fig
% 
% subplot(1,2,2)
% boxplot(100*Miss_ADHD, 'Colors',Colors(2,:), 'symbol', 'o k')
% ylim([0 15])
% ylabel('Misses')
% xlabel('ADHD')
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     set(h(j),'LineWidth',lineWidth)
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.5);
% end
% format_fig
% 
% 
% % FAs
% lineWidth = 3;
% figure; 
% subplot(1,2,1)
% boxplot(100*FA_CTR,'Colors',Colors(1,:), 'symbol', 'o k')
% ylim([0 100])
% ylabel('False Alarms')
% xlabel ('Control')
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     set(h(j),'LineWidth',lineWidth)
% 
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.5);
% end
% format_fig
% 
% subplot(1,2,2)
% boxplot(100*FA_ADHD, 'Colors',Colors(2,:), 'symbol', 'o k')
% ylim([0 100])
% ylabel('False Alarms')
% xlabel('ADHD')
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     set(h(j),'LineWidth',lineWidth)
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.5);
% end
% format_fig


%% Mind States 
numbers_of_interest = [1, 2, 3, 4];
ADHD_state = [];
Control_state =[];
Control_state = zeros(length(ctrs), length(numbers_of_interest));
ADHD_state = zeros(length(adhds), length(numbers_of_interest));

clear state_values
for nc = 1:length(ctrs)
    % Extract the State column for the current participant
    state_values = all_probe_table.State(all_probe_table.SubID == ctrs(nc));
    % Count occurrences of each number (1, 2, 3, 4) for the current participant
    Control_state(nc, :) = histcounts(state_values, [numbers_of_interest, Inf]);
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
labels = {'On Task', 'Mind Wandering','Mind Blanking', 'Dont Remember'};

figure('Position',[1955         361         758         413]);
bar((1:numel(labels))-0.2, Ctr_state_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',0.38);
xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);
% Add percentage values on top of each bar
for i = 1:numel(labels)
    for j = 1:size(Ctr_state_percentage_distribution, 1)
        text(i-0.2, Ctr_state_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_state_percentage_distribution(j, i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12,'Color',Colors(1,:),'FontWeight','bold');
    end
end
format_fig
hold on;
bar((1:numel(labels))+0.2, ADHD_state_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',0.38);
ylabel('% of Mind State');
xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);
% Add percentage values on top of each bar
for i = 1:numel(labels)
    for j = 1:size(ADHD_state_percentage_distribution, 1)
        text(i+0.2, ADHD_state_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_state_percentage_distribution(j, i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12,'Color',Colors(2,:),'FontWeight','bold');
    end
end
format_fig
xlim([0.5 4.5])

%% Vigilance
numbers_of_interest = [1, 2, 3, 4];
ADHD_vig = [];
Control_vig =[];
Control_vig = zeros(length(ctrs), length(numbers_of_interest));
ADHD_vig = zeros(length(adhds), length(numbers_of_interest));

clear vig_values
for nc = 1:length(ctrs)
    % Extract the State column for the current participant
    vig_values = all_probe_table.Vigilance(all_probe_table.SubID == ctrs(nc));
    % Count occurrences of each number (1, 2, 3, 4) for the current participant
    Control_vig(nc, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
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



% Plot the distribution
labels = {'Extremely Alert', 'Alert','Sleepy', 'Extremely Sleepy'};

figure('Position',[1955         361         758         413]);
bar((1:numel(labels))-0.2, Ctr_vig_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',0.38);
ylabel('% Vigilance Ratings');
xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);
% Add percentage values on top of each bar
for i = 1:numel(labels)
    for j = 1:size(Ctr_vig_percentage_distribution, 1)
        text(i-0.2, Ctr_vig_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_vig_percentage_distribution(j, i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 12,'Color',Colors(1,:),'FontWeight','bold');
    end
end
format_fig
hold on;

bar((1:numel(labels))+0.2, ADHD_vig_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',0.38);
xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);
% Add percentage values on top of each bar
for i = 1:numel(labels)
    for j = 1:size(ADHD_vig_percentage_distribution, 1)
        text(i+0.2, ADHD_vig_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_vig_percentage_distribution(j, i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12,'Color',Colors(2,:),'FontWeight','bold');
    end
end
format_fig
xlim([0.5 4.5])

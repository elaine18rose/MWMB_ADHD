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

%% Loading data 


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
    block_table=zeros(4,15);
    block_table=array2table(block_table,...
        'VariableNames',{'SubID','Group','Block','ON','MW','MB','DK','Dist','Perso','Int','Intention','Vigilance','Miss','FA','HitRT'});
    block_table.SubID=categorical(block_table.SubID);
    block_table.Group=categorical(block_table.Group);

    block_table.BlockN=(1:4)';
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




%%
%%%%% EXAMPLE OF STATS
% trial-level performance
all_behav_table.Group=categorical(all_behav_table.Group);
all_behav_table.SubID=categorical(all_behav_table.SubID);
mdl=fitlme(all_behav_table,'FA~1+BlockN*Group+(1|SubID)'); % you can do this for miss and RT

% block-level stats
all_block_table.Group=categorical(all_block_table.Group);
all_block_table.SubID=categorical(all_block_table.SubID);
mdl=fitlme(all_block_table,'FA~1+BlockN*Group+(1|SubID)');
mdl=fitlme(all_block_table,'ON~1+BlockN*Group+(1|SubID)');
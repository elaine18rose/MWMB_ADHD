%% ERP for trials - changes from v1: re-referencing and filtering again

clear all;
close all;

if isempty(findstr(pwd,'thandrillon'))==0
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    data_path='/Users/thandrillon/Data/ADHD_MW/EEG/';
    behav_path = '/Users/thandrillon/Data/ADHD_MW/Behaviour/';
    preproc_path='/Users/thandrillon/Data/ADHD_MW/Preproc/';
    path_eeglab='/Users/thandrillon/WorkGit/projects/ext/eeglab/';
    path_ICAlabel='/Users/thandrillon/WorkGit/projects/ext/ICLabel/';
    path_FMINSEARCHBND='/Users/thandrillon/Work/local/FMINSEARCHBND';
    path_ExGauss='/Users/thandrillon/WorkGit/projects/ext/exgauss/';
else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/desktop/MATLAB_Functions/fieldtrip/';
    data_path = '/Volumes/Seagate/MWMB_ADHD_SART/EEG/';
    behav_path = '/Volumes/Seagate/MWMB_ADHD_SART/Behaviour/';
    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    path_detectSW = '/Volumes/Seagate/MWMB_ADHD_SART/SW_detection/';
    path_eeglab='/Users/elaine/desktop/MATLAB_Functions/eeglab/';
    path_ICAlabel='/Users/elaine/desktop/MATLAB_Functions/ICLabel/';
    path_ExGauss='/Users/elaine/desktop/MATLAB_Functions/exgauss';
    path_FMINSEARCHBND='/Users/elaine/desktop/MATLAB_Functions/FMINSEARCHBND';
    path_RainCloudPlot='/Users/elaine/desktop/MATLAB_Functions/RainCloudPlots/';

    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(path_LSCPtools))
% addpath(genpath(path_RainCloudPlot));
addpath(path_fieldtrip)
ft_defaults;
addpath(genpath(path_ExGauss))
addpath(genpath(path_FMINSEARCHBND))



% NOTE: I'm using trial data 
files=dir([preproc_path filesep 'clean_i_trial_*.mat']); % EP - changed 14/11/24 files=dir([preproc_path filesep 'fetrial_ft_*.mat']);

load([pwd filesep '..' filesep 'Preproc' filesep 'all_badChannels_badProbes.mat']);

%EEG Layout info
run ../MWMB_ADHD_elec_layout.m

%% Loop across files
redo=0;

nFc=0;
nFc4=0;
all_ERP_NG=[];
all_ERP_G=[];

group_PowDataEO=[];

for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(15:end-4); 

%     if ismember(SubID, {'A053', 'C012', 'C024', 'C036', 'C038'}) %% TEMPORARY SKIP CAUSE IT WAS CRASHING; C012 & C024 & C038 maybe just rerun importSART
%         continue;
%     end

    behav_files = dir([behav_path filesep 'wanderIM_behavres_' SubID '_*.mat']);
    if isempty(behav_files)
        warning('Cannot find the behavioral file for %s\n', SubID);
        continue;
    end
    if exist([preproc_path  filesep 'ERP_' file_name])==0 || redo==1 || nF==1 
        table=load([behav_path behav_files.name]);  
        orifile=dir([data_path filesep 'ID-' SubID '.eeg']);
        
        nFc=nFc+1;
        fprintf('... processing %s\n',file_name);
        load([preproc_path file_name]);

%         %%% rename channels
%         mylabels = data.label;
%         for nCh = 1:length(mylabels)
%             findspace = findstr(mylabels{nCh},' ');
%             if isempty(findspace)
%                 newlabels{nCh} = mylabels{nCh};
%             else
%                 if ismember(mylabels{nCh}(1),{'1','2','3','4','5','6','7','8','9'})
%                     newlabels{nCh} = mylabels{nCh}(findspace+1:end);
%                 else
%                     newlabels{nCh} = mylabels{nCh}(1:findspace-1);
%                 end
%             end
%         end
%         cfg=[];
%         cfg.channel          = data.label; %hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG'))) & (cellfun(@isempty,regexp(hdr.label,'EMG'))) & (cellfun(@isempty,regexp(hdr.label,'ECG'))) & (cellfun(@isempty,regexp(hdr.label,'Mic')))));
%         cfg.montage.labelold = data.label;
%         cfg.montage.labelnew = newlabels;
%         cfg.montage.tra      = eye(length(data.label));
%         data = ft_preprocessing(cfg, data);
       
        cfg=[];   %% !! EP - put back in from v1 although should this go after baseline correcting?
        cfg.reref           = 'yes'; 
        cfg.refchannel      = {'TP9', 'TP10'}; %'all';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        % cfg.dftfreq        = [25]; % set up the frequencies for notch filtering; EP - left this commented out since there's no flicker
        data_clean = ft_preprocessing(cfg,data);
        
        
        %%% take out trial %% !!! EP - I commented this out as we took out probes not trials 
%         thisF=match_str(badChannels_table(:,1),SubID); 
%         if isempty(thisF)
%             warning('cannot find the correct line in the bad channels table\n')
%             continue;
%         end
%         badTrials = badChannels_table{thisF,2}; %eval(['badTrials=[' badChannels_table.Bad_Trials{thisF} '];']);
%         old_trials=old_table.Sample(badTrials);
%         table.StimType(ismember(table.Sample,old_trials))=NaN;
%         table.StimType(find(diff(table.BlockN)==1))=NaN;
        
        % COMPUTE BOTH ONSET AND OFFSET ERP
        cfgerp        = [];
        %     cfgerp.latency        = [-0.2 0.8];
        cfgerp.trials = find_trials({events.value},'S 10'); % EP  - Here looking for NoGos
        av_data_NG = ft_timelockanalysis(cfgerp, data_clean); % av_data_NG = ft_timelockanalysis(cfgerp, data);

        cfgerp        = [];
        cfgerp.trials = find_trials({events.value},'S  9'); % EP - Here looking for Gos
        av_data_G = ft_timelockanalysis(cfgerp, data_clean); % av_data_G = ft_timelockanalysis(cfgerp, data);

        ERP_NG=av_data_NG.avg-repmat(nanmean(av_data_NG.avg(:,av_data_NG.time>-0.2 & av_data_NG.time<0),2),1,length(av_data_NG.time));
        % ERP_G=av_data_G.avg-repmat(nanmean(av_data_G.avg(:,av_data_NG.time>-0.2 & av_data_NG.time<0),2),1,length(av_data_NG.time)); % This is using av_data_NG for baseline correcting ERP_G
        ERP_G=av_data_G.avg-repmat(nanmean(av_data_G.avg(:,av_data_G.time>-0.2 & av_data_G.time<0),2),1,length(av_data_G.time)); % TA - This is doing essentially the same as the above 
       
        save([preproc_path  filesep 'ERP_' file_name(1:end-4)],'ERP_NG','ERP_G')
    else
        load([preproc_path filesep 'ERP_' file_name])
          nFc=nFc+1;
    end
    
    all_ERP_NG(nFc,:,:)=ERP_NG;
    all_ERP_G(nFc,:,:)=ERP_G;
  
    
    if SubID(1)=='C'
        group_PowDataEO{nFc}='Control';
    elseif SubID(1)=='A'
        group_PowDataEO{nFc}='ADHD';
    end

    SubIDs{nFc} = SubID;

end

%%
xTime=av_data_NG.time;
chLabels=av_data_NG.label;

diff_all_ERP=all_ERP_G-all_ERP_NG;

%% Plots
%Event-related potentials non-target vs target

Colors1=[233,163,201;
    247,247,247;
    161,215,106]/256;

figure;
thisCh=match_str(chLabels,'Fz');
subplot(2,2,1);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(:,thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(:,thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'No Go','Go'})
title('ERP NoGo vs Go: Fz');
xlim([-0.2 1.6])
format_fig;

thisCh=match_str(chLabels,'Cz');
subplot(2,2,2);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(:,thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(:,thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'No Go','Go'})
title('ERP NoGo vs Go: Cz');
xlim([-0.2 1.6])
format_fig;

thisCh=match_str(chLabels,'Pz');
subplot(2,2,3);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(:,thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(:,thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'No Go','Go'})
title('ERP NoGo vs Go: Pz');
xlim([-0.2 1.6])
format_fig;

thisCh=match_str(chLabels,'Oz');
subplot(2,2,4);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(:,thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(:,thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'No Go','Go'})
title('ERP NoGo vs Go: Oz');
xlim([-0.2 1.6])
format_fig;


%% ERP and split by group - manuscript figure
f1=figure('Position', [100, 100, 1000, 1000]); 
thisCh=match_str(chLabels,'Fz');
subplot(2,2,1);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(3,:),0,':',0.5,1,0,1,2);
hold on;
%legend(hp,{'No Go','Go'})
title('Fz');
xlabel('Time (s)')
ylabel('Voltage (µV)')
ylim([-4 8])
xlim([-0.2 1.6])
format_fig;
disp(['Min ERP value: ', num2str(min(all_ERP_NG(:)))]);
disp(['Max ERP value: ', num2str(max(all_ERP_NG(:)))]);

thisCh=match_str(chLabels,'Cz');
subplot(2,2,2);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(3,:),0,':',0.5,1,0,1,2);
hold on;
lgd = legend(hp,{'NoGo Control','Go Control','NoGo ADHD','Go ADHD'},'Location','northeast', 'Box', 'off');
lgd.Position = [0.85, 0.8, 0.1, 0.1];
title('Cz');
xlabel('Time (s)')
ylabel('Voltage (µV)')
ylim([-4 15])
xlim([-0.2 1.6])
format_fig;

thisCh=match_str(chLabels,'Pz');
subplot(2,2,3);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(3,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'No Go','Go'})
xlabel('Time (s)')
ylabel('Voltage (µV)')
ylim([-4 15])
title('Pz');
xlim([-0.2 1.6])
format_fig;

thisCh=match_str(chLabels,'Oz');
subplot(2,2,4);
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(3,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'No Go','Go'})
xlabel('Time (s)')
ylabel('Voltage (µV)')
ylim([-4 8])
title('Oz');
xlim([-0.2 1.6])
format_fig;

sgtitle('ERP NoGo vs Go', 'FontWeight', 'bold', 'FontSize', 25);
% Save figure
% saveas(gcf, [pwd filesep 'Figures' filesep 'Fig3_PanelA_ERP.svg']); 

%% ERP differences
figure;
thisCh=match_str(chLabels,'Fz'); 
subplot(2,2,1); 
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'Non Target','Target'})
title('ERP differences NoGo vs Go split by group: Fz');
xlim([-0.2 1.5])
format_fig;

thisCh=match_str(chLabels,'Cz'); 
subplot(2,2,2); 
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'Non Target','Target'})
title('ERP differences NoGo vs Go split by group: Cz');
xlim([-0.2 1.5])
format_fig;

thisCh=match_str(chLabels,'Pz'); 
subplot(2,2,3); 
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'Non Target','Target'})
title('ERP differences NoGo vs Go split by group: Pz');
xlim([-0.2 1.5])
format_fig;

thisCh=match_str(chLabels,'Oz'); 
subplot(2,2,4); 
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors1(1,:),0,':',0.5,1,0,1,2);
hold on;
% legend(hp,{'Non Target','Target'})
title('ERP differences NoGo vs Go split by group: Oz');
xlim([-0.2 1.5])
format_fig;


%% Stats for Manuscript part 1:
thisCh=match_str(chLabels,'Pz'); %Pz, Fz
Group_A=squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:));
Group_B=squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:));
Groups=[ones(size(Group_A,1),1) ; 2*ones(size(Group_B,1),1)];
totPerm=1000;
[realpos]=get_cluster_permutation_aov([Group_A ; Group_B],Groups,0.1,0.05,totPerm,xTime,[],[]); % Runs an ANOVA - Main effect of group on the diff


%% Stats for Manuscript part 2:
Colors2 = [175, 175, 175;  % light grey (darker than off-white)
    247, 247, 247;  % off-white
    161, 215, 106   % light green
] / 256;
f1=figure('Position', [100, 100, 1000, 1000]); 
thisCh=match_str(chLabels,'Pz');
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors2(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors2(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:)),0,Colors2(1,:),0,':',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(xTime,squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,Colors2(3,:),0,':',0.5,1,0,1,2);
hold on;
lgd = legend(hp,{'NoGo Control','NoGo ADHD','Go Control','Go ADHD'},'Location','northeast', 'Box', 'off');

Group_A=squeeze(all_ERP_NG(match_str(group_PowDataEO,'Control'),thisCh,:))-squeeze(all_ERP_G(match_str(group_PowDataEO,'Control'),thisCh,:));
Group_B=squeeze(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),thisCh,:))-squeeze(all_ERP_G(match_str(group_PowDataEO,'ADHD'),thisCh,:));

both_Groups=[Group_A ; Group_B];
for k=1:size(both_Groups,2)
    [p_aov(k),anovatab,stats] =anova1(both_Groups(:,k),Groups,'off');
    [p_ttest2(k)] =ranksum(Group_A(:,k),Group_B(:,k));
end
% line(xTime(find(p_aov<0.05)),14*ones(1,sum(p_aov<0.05)),'Color','k','LineWidth',2,'LineStyle',':') % code from TA

ylimits = ylim; % Get y-axis limits
sig_time = find(p_aov < 0.05); % Find significant frequencies

% shading area that is sig. different
if ~isempty(sig_time)
    sig_time_values = xTime(sig_time); % Get corresponding frequency values
    hp(5) = fill([sig_time_values, fliplr(sig_time_values)], ...
         [ylimits(1)*ones(1, length(sig_time_values)), ylimits(2)*ones(1, length(sig_time_values))], ...
         [0.65 0.65 0.65], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Significant Diff');
end

if ~isempty(sig_time_values)
    sig_range = [min(sig_time_values), max(sig_time_values)];
    disp(['Significant time range: ', num2str(sig_range(1)), ' - ', num2str(sig_range(2)), ' s']);
else
    disp('No significant time range found.');
end

legend;
xlabel('Time (s)')
ylabel('Voltage (µV)')
ylim([-4 15])
xlim([-0.2 1.6])
format_fig;
lgd = legend(hp,{'NoGo NT','NoGo ADHD','Go NT','Go ADHD', 'Significant Diff'},'Location','northeast', 'Box', 'off', 'FontSize', 32);
lgd.Position = [0.67, 0.72, 0.1, 0.1]; %[X, Y, Width, Height]
set(gca, 'FontSize', 35);
title('ERP NoGo vs Go - Pz', 'FontWeight', 'bold', 'FontSize', 45);

saveas(gcf, [pwd filesep 'Figures' filesep 'Fig3_PanelA_ERP.svg']); 



%% Stats 
TimeWin=[0.15 0.25];
thisCh=match_str(chLabels,'POz');
diffG_NG=squeeze(mean(all_ERP_G(:,thisCh,xTime>TimeWin(1) & xTime<TimeWin(2))-all_ERP_NG(:,thisCh,xTime>TimeWin(1) & xTime<TimeWin(2)),3));

Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(diffG_NG(match_str(group_PowDataEO,'Control')), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',1,'bound_data',[0 100]);
h2 = raincloud_plot(diffG_NG(match_str(group_PowDataEO,'ADHD')), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',1,'bound_data',[0 100]);
set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
%bset(gca,'XLim', [-30 40], 'YLim', ylim.*[0 0.05]);
format_fig; title('Diff Go/NoGos (s)'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});


[h, pV_diffGroup,~,stat_diffGroup]=ttest2(diffG_NG(match_str(group_PowDataEO,'Control')),diffG_NG(match_str(group_PowDataEO,'ADHD')));
fprintf('... unpaired t-test between groups on diff ERP between Go and NoGo on POz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) 

[h, pV_CTR,~,stat_CTR]=ttest(diffG_NG(match_str(group_PowDataEO,'Control')),0);
fprintf('... one-sample t-test with 0 for Controls on diff ERP between Go and NoGo on POz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_CTR,stat_CTR.tstat,stat_CTR.df)
[h, pV_ADHD,~,stat_ADHD]=ttest(diffG_NG(match_str(group_PowDataEO,'ADHD')),0);
fprintf('... unpaired t-test with 0 for ADHD on diff ERP between Go and NoGo on POz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_ADHD,stat_ADHD.tstat,stat_ADHD.df)

% %%
% %Event-related potentials non-target vs target OFFSET
% thisCh=match_str(chLabels,'Fz');
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NT_offset(:,thisCh,:)),0,'b',0,'-',0.5,1,0,1,2);
% hold on;
% [~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_TG_offset(:,thisCh,:)),0,'r',0,'-',0.5,1,0,1,2);
% hold on;
% legend(hp,{'Non Target','Target'})
% title('Event-related potentials non-target vs target OFFSET');
% format_fig;

%%
%%% Plot topography [0.1-0.3]s post-offset
% cfg = [];
% cfg.layout = 'eeg1010.lay';
% layout=ft_prepare_layout(cfg);

cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.channel = data.label;
cfg.channel=data.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

matching_elec=[];
for nE=1:length(layout.label)-2
    matching_elec(nE)=(match_str(data.label,layout.label(nE)));
end

%Difference Go/NoGo trials onset-locked
temp_topo=squeeze(mean(mean(diff_all_ERP(:,matching_elec,xTime>TimeWin(1) & xTime<TimeWin(2)),3),1));
figure;
subplot(2,2,1)
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Topography difference Go/NoGo trials [1.1-1.8]s post-onset')
caxis([-1 1]*2.5)
format_fig;

temp_topo_CTR=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.15 & xTime<0.25),3),1));
subplot(2,2,2)
simpleTopoPlot_ft(temp_topo_CTR', layout,'on',[],0,1);
colorbar;
title('Topography difference Go/NoGo trials [1.1-1.8]s post-onset for Controls')
caxis([-1 1]*2.5)
format_fig;

temp_topo_ADHD=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.15 & xTime<0.25),3),1));
subplot(2,2,3)
simpleTopoPlot_ft(temp_topo_ADHD', layout,'on',[],0,1);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset for ADHDs')
caxis([-1 1]*2.5)
format_fig;

%T-values topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(diff_all_ERP,2)
    A=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),nE,xTime>0.15 & xTime<0.25),3));
    B=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),nE,xTime>0.15 & xTime<0.25),3));
    [h,pV,~,stat]=ttest2(B,A);
    temp_topo_tval(nE)=stat.tstat;
    temp_topo_pval(nE)=pV;
end
subplot(2,2,4)
simpleTopoPlot_ft(temp_topo_tval(matching_elec), layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset (tvalue)')
caxis([-1 1]*2.5)
if ~isempty(find(temp_topo_pval(matching_elec)<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval(matching_elec)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
format_fig;
%% below I think we can change it
%%% Plot topography [0.3-0.8]s post-offset
% %Non-target trials offset-locked
% temp_topo=squeeze(mean(mean(all_ERP_NG_offset(:,matching_elec,xTime>TimeWin(1) & xTime<TimeWin(2)),3),1));
% figure;
% simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
% colorbar;
% title('Topography NoGo trials [0.3-0.8]s post-offset')
% caxis([-1 1])
% %Target trials offset-locked
% temp_topo=squeeze(mean(mean(all_ERP_G_offset(:,matching_elec,xTime>TimeWin(1) & xTime<TimeWin(2)),3),1));
% figure;
% simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
% colorbar;
% title('Topography target trials [0.3-0.8]s post-offset')
% caxis([-1 1])
% %Difference TG/NT trials offset-locked
% temp_topo=squeeze(mean(mean(diff_all_ERP_offset(:,matching_elec,xTime>TimeWin(1) & xTime<TimeWin(2)),3),1));
% figure;
% simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
% colorbar;
% title('Topography difference target/non target trials [0.3-0.8]s post-offset')
% caxis([-1 1])

%% Difference TG/NG for all subjects

thisCh=match_str(chLabels,'Oz');
figure;
[~,~]=simpleTplot(xTime,squeeze(diff_all_ERP(:,thisCh,:)),0,'k',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
title('ERP differences between Target and Non-target trials for all subjects ONSET');

%%
%Difference TG/NG for Controls and ADHDs
thisCh=match_str(chLabels,'Fz');

figure;
hp=[];
% [~,hp(1)]=simpleTplot(xTime,squeeze(diff_all_ERP(:,thisCh,:)),0,'k');
% hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(diff_all_ERP(match_str(group_PowDataEO,'Control'),thisCh,:)),0,'b',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,'r',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('ERP differences between Target and Non-target trials ONSET');


Group_A=squeeze(diff_all_ERP(match_str(group_PowDataEO,'Control'),thisCh,:));
Group_B=squeeze(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),thisCh,:));
Groups=[ones(size(Group_A,1),1) ; 2*ones(size(Group_B,1),1)];
%[realpos realneg]=get_cluster_permutation_aov([Group_A ; Group_B],Groups,0.2,0.1,1000,xTime_offset,'full',[]); % Runs an ANOVA - Main effect of group on the diff


%% Group topographies [0.2-0.25]s
% NoGo trials - Controls
temp_topo=squeeze(mean(mean(all_ERP_NG(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
figure;
subplot(2,3,1)
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-1 5])
title('Controls NoGo 0.2-0.25s')
format_fig
hold on

% Go trials - Controls
temp_topo=squeeze(mean(mean(all_ERP_G(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
subplot(2,3,2)
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-1 5])
title('Controls Go 0.2-0.25s')
format_fig
hold on

%Difference G/NG trials - Controls
temp_topo=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
subplot(2,3,3);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-3 3])
title('Controls - diff Go/NoGo 0.2-0.25s')
format_fig
hold on

% NoGo trials - ADHDs
temp_topo=squeeze(mean(mean(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
subplot(2,3,4);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-1 5])
title('ADHD NoGo 0.2-0.25s')
format_fig
hold on

% Go trials - ADHDS
temp_topo=squeeze(mean(mean(all_ERP_G(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
subplot(2,3,5);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-1 5])
title('ADHD Go 0.2-0.25s')
format_fig
hold on

% Difference Go/NoGo trials - ADHDs
temp_topo=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.2 & xTime<0.25),3),1));
subplot(2,3,6);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([-3 3])
title('ADHD - diff Go/NoGo 0.2-0.25s')
format_fig
hold off

%% Group topographies [0.4-0.5]s
% NoGo trials - Controls
temp_topo=squeeze(mean(mean(all_ERP_NG(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
figure;
subplot(2,3,1)
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([0 10])
title('Controls NoGo 0.4-0.5s')
format_fig
hold on

% Go trials - Controls
temp_topo=squeeze(mean(mean(all_ERP_G(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
subplot(2,3,2)
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([0 4])
title('Controls Go 0.4-0.5s')
format_fig
hold on

%Difference G/NG trials - Controls
temp_topo=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
subplot(2,3,3);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
%caxis([-5 5])
title('Controls - diff Go/NoGo 0.4-0.5s')
format_fig
hold on

% NoGo trials - ADHDs
temp_topo=squeeze(mean(mean(all_ERP_NG(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
subplot(2,3,4);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([0 10])
title('ADHD NoGo 0.4-0.5s')
format_fig
hold on

% Go trials - ADHDS
temp_topo=squeeze(mean(mean(all_ERP_G(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
subplot(2,3,5);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
caxis([0 4])
title('ADHD Go 0.4-0.5s')
format_fig
hold on

%Difference Go/NoGo trials - ADHDs
temp_topo=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.4 & xTime<0.5),3),1));
subplot(2,3,6);
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
%caxis([-5 5])
title('ADHD - diff Go/NoGo 0.4-0.5s')
format_fig
hold off

%% tvalue calc and plots for 0.2-0.25s and 0.4-0.5s
temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(diff_all_ERP,2)
    A=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),nE,xTime>0.2 & xTime<0.25),3));
    B=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),nE,xTime>0.2 & xTime<0.25),3));
    [h,pV,~,stat]=ttest2(B,A);
    temp_topo_tval(nE)=stat.tstat;
    temp_topo_pval(nE)=pV;
end
figure
subplot(1,2,1)
simpleTopoPlot_ft(temp_topo_tval(matching_elec), layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Difference Go/NoGo 0.2-0.25s (tvalue)')
caxis([-1 1]*3)
if ~isempty(find(temp_topo_pval(matching_elec)<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval(matching_elec)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
format_fig;
hold on 

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(diff_all_ERP,2)
    A=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),nE,xTime>0.4 & xTime<0.5),3));
    B=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),nE,xTime>0.4 & xTime<0.5),3));
    [h,pV,~,stat]=ttest2(B,A);
    temp_topo_tval(nE)=stat.tstat;
    temp_topo_pval(nE)=pV;
end
subplot(1,2,2)
simpleTopoPlot_ft(temp_topo_tval(matching_elec), layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Difference Go/NoGo 0.4-0.5s (tvalue)')
caxis([-1 1]*3)
if ~isempty(find(temp_topo_pval(matching_elec)<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval(matching_elec)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
format_fig;
hold off 

%% Cluster difference for ADHD and controls 
temp_topoADHD=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>0.1 & xTime<0.3),3),1));
temp_topoCTRL=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),matching_elec,xTime>0.1 & xTime<0.3),3),1));

temp_topo = temp_topoADHD-temp_topoCTRL; %usually good to do group of interest minus control condition
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Cluster difference for ADHD and controls [0.1-0.3]s')


%%
thisCh=match_str(chLabels,'Fz');

All_Conds=double(ismember(group_PowDataEO,'Control'))+1;
% [realpos_lin realneg_lin]=get_cluster_permutation_aov(squeeze(diff_all_ERP(:,thisCh,:)),All_Conds',...
%     0.2,0.1,1000,xTime); % Commented this out because it's taking too
%     long to run



%% Checking if there's a sig difference between P300 for NoGo and Go 
G_data = squeeze(all_ERP_G(:,25,:)); % 25 because it's Pz 
G_data_p300= G_data(:,125); % 125 because it's where P300 occurs; to do this, Mana looked at the sampling frequency and examined the plot I showed her (how it went from -2 to 800 ms) and input it in a script
%NOTE: The note in the previous line is for CTET data. May need to modify this for SART data
NG_data = squeeze(all_ERP_NG(:,25,:));
NG_data_p300= NG_data(:,125);
[h p] = ttest(NG_data_p300, G_data_p300) % T test to see if there's a difference between NT and TG







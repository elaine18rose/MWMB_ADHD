%% version 2 - uses version 2 getSW table with NA values
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

    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(path_LSCPtools))
addpath(path_fieldtrip)
ft_defaults;
addpath(genpath(path_ExGauss))
addpath(genpath(path_FMINSEARCHBND))

% select relevant files, here task
SW_files=dir([preproc_path filesep 'SW_clean_i_probe_*.mat']);

%EEG Layout info
run ../MWMB_ADHD_elec_layout.m


%% Loop across files
% cfg = [];
% cfg.layout = 'EEG1010.lay';
% cfg.channel=cellstr(Elec_group);
% cfg.center      = 'yes';
% layout=ft_prepare_layout(cfg);


%load([pwd filesep '..' filesep 'Preproc' filesep 'all_badChannels_badProbes.mat']);
load([pwd filesep '..' filesep 'Preproc' filesep 'all_badChannels_badProbes.mat']);

%%
SW_table=array2table(zeros(0,17),'VariableNames',{'SubID','Group','Block','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold','SW_peakneg','SW_peakpos','Probe_MS','Probe_Vig','Behav_Miss','Behav_FA','Behav_RT'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Group=categorical(SW_table.Group);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.Probe_MS=categorical(SW_table.Probe_MS);

MS_labels={'ON','MW','MB','DK'};
FilesPbme=[];
bins=0:2:150;
all_distrib_SW=[];
all_distrib_SW_per=[];
all_GroupID=[];
for nF=1:length(SW_files)
    if startsWith(SW_files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end


    %%% load the data
    SubInfo=split(SW_files(nF).name,'_');
    SubID=SubInfo{end}(1:end-4);
    if SubID(1)=='A'
        GroupID='ADHD';
    elseif SubID(1)=='C'
        GroupID='Control';
    else
        GroupID='undefined';
    end

    if strcmp(SubID,'C017')
        FilesPbme=[FilesPbme ; {SubID} , {'Noisy SW detection'} , {[]}];
        continue;
    end
    if strcmp(SubID,'C015')
        FilesPbme=[FilesPbme ; {SubID} , {'Weird behaviour'} , {[]}];
        continue;
    end

    behav_file=dir([data_path filesep '..' filesep 'Behaviour' filesep 'wanderIM_behavres_' SubID '*.mat']);
    if length(behav_file)
        load([behav_file.folder filesep behav_file.name])
    else
        FilesPbme=[FilesPbme ; {SubID} , {'Missing Behaviour'}, {[]}];
        continue;
    end

    fsample=500;
    labels={'Fp1'	'Fp2'	'F7'	'F3'	'Fz'	'F4'	'F8'	'FC5'	'FC1'	'FC2'	'FC6'	'T7'	'C3'	'Cz'	'C4'	'T8'	'TP9'	'CP5'	'CP1'	'CP2'	'CP6'	'TP10'	'P7'	'P3'	'Pz'	'P4'	'P8'	'PO9'	'O1'	'Oz'	'O2'	'PO10'	'AF7'	'AF3'	'AF4'	'AF8'	'F5'	'F1'	'F2'	'F6'	'FT9'	'FT7'	'FC3'	'FC4'	'FT8'	'FT10'	'C5'	'C1'	'C2'	'C6'	'TP7'	'CP3'	'CPz'	'CP4'	'TP8'	'P5'	'P1'	'P2'	'P6'	'PO7'	'PO3'	'POz'	'PO4'	'PO8'};
    load([preproc_path filesep 'SW_clean_i_probe_' SubID]);

    if ~isempty(match_str(badChannels_badTrials_info(:,1),SubID)) && ~isempty(badChannels_badTrials_info{match_str(badChannels_badTrials_info(:,1),SubID),7})
        probe_res(badChannels_badTrials_info{match_str(badChannels_badTrials_info(:,1),SubID),7},:)=[];
    elseif strcmp(SubID,'A008')
        probe_res(25,:)=[];
    elseif strcmp(SubID,'A053')
        probe_res([24 25],:)=[]; %removing two probes because we cannot find the triggers
    elseif strcmp(SubID,'C036')
        probe_res([6 7],:)=[]; %removing two probes because we cannot find the triggers
    end
    %     if size(probe_res,1)~=length(unique(all_Waves(:,2)))
    %         FilesPbme=[FilesPbme ; {SubID} , {'Different Numbers of Probes'}, {[size(probe_res,1) length(unique(all_Waves(:,2)))]}];
    %         continue;
    %     end

    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[4 7]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;

    all_Waves=double(all_Waves);
    all_onset=all_Waves(:,5)/fsample-25;
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./fsample);
    fprintf('... ... %g %% waves discarded because of onset\n',mean(all_onset<-20 | all_onset>0)*100)
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_onset<-20 | all_onset>0 | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];


     for nBl=unique(all_Waves(:,2))' % EP - uncommmented
        duration_of_recording=20/60*length(unique(all_Waves(:,2)));

        table_length=size(SW_table,1);
        SW_table.SubID(table_length+1)={SubID};
        SW_table.Group(table_length+1)={GroupID};
        SW_table.Block(table_length+1)=nBl; % EP - uncommmented 141 - 156


        this_MS=MS_labels(probe_res(nBl,19));
        SW_table.Probe_MS(table_length+1)=this_MS;
        SW_table.Probe_Vig(table_length+1)=probe_res(nBl,22);

        this_block=probe_res(nBl,4);
        this_trial=probe_res(nBl,6);
        temp_gos=test_res(test_res(:,1)==this_block & test_res(:,4)<this_trial & test_res(:,5)~=3,:);
        temp_gos=temp_gos(end-17:end,:);
        temp_nogos=test_res(test_res(:,1)==this_block & test_res(:,4)<this_trial & test_res(:,5)==3,:);
        temp_nogos=temp_nogos(end-1:end,:);
        SW_table.Behav_Miss(table_length+1)=nanmean(temp_gos(:,12)==0);
        SW_table.Behav_RT(table_length+1)=nanmean(temp_gos(:,10)-temp_gos(:,8));
        SW_table.Behav_FA(table_length+1)=nanmean(temp_nogos(:,11)==0);
     end % EP - uncommmented

        this_amp=all_Waves(:,4);
        amp_bins=histc(this_amp,bins);
        amp_bins=(amp_bins/duration_of_recording)/length(labels)';

        all_distrib_SW=[all_distrib_SW ; amp_bins'];
        temp_distr_SW=[];
        for nE=1:length(labels)
            this_amp=all_Waves(all_Waves(:,3)==nE,4);
            amp_bins=histc(this_amp,bins);
            amp_bins=(amp_bins/duration_of_recording)';
            temp_distr_SW=[temp_distr_SW ; amp_bins];
        end
                all_distrib_SW_per=cat(3,all_distrib_SW_per, temp_distr_SW);
        all_GroupID=[all_GroupID; {GroupID}];


end


%%
figure('Position',[347   308   441   306]);
uniqueGroupIDs={'Control','ADHD'};
uniqueGroupLabel={'Neurotypical','ADHD'};
all_GroupID=categorical(all_GroupID);
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;
hplot=[];
for nG=1:length(uniqueGroupIDs)
    temp_plot=all_distrib_SW(all_GroupID==uniqueGroupIDs(nG),:);
    [pV hplot(nG)]=simpleTplot(bins,temp_plot,0,Colors(nG,:),[0],'-',0.5,1,0,1,4);
end
legend(hplot,uniqueGroupLabel, 'Box', 'off')
format_fig;
xlabel('Amplitude (µV)')
ylabel('Density (wave/min)')
xlim([0 75])
% Save figure
saveas(gcf, [pwd filesep 'Figures' filesep 'Fig4_PanelAi_Dist.svg']);

figure('Position',[347   308   441   306]);
uniqueGroupIDs={'Control','ADHD'};
hplot=[];
% for nG=1:length(uniqueGroupIDs)
    temp_plotA=all_distrib_SW(all_GroupID==uniqueGroupIDs(2),:);
     temp_plotB=all_distrib_SW(all_GroupID==uniqueGroupIDs(1),:);
     temp_plotA=temp_plotA-nanmean(temp_plotB);
   simpleTplot(bins,temp_plotA,0,Colors(2,:),[0 0.05 0.05 1000],'-',0.5,1,0,1,4);
% end
% legend(hplot,uniqueGroupIDs)
yline(0, 'k', 'LineWidth', 2);
format_fig;
xlabel('Amplitude (µV)')
ylabel({'Density Diff (ADHD - NT)', 'waves/min'});
xlim([0 75])
ylim([-1 1]*3)
% Save figure
saveas(gcf, [pwd filesep 'Figures' filesep 'Fig4_PanelAii_Dist.svg']);

%% Per electrode
figure('Position',[347   308   441   306]);
uniqueGroupIDs={'Control','ADHD'};
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;
hplot=[];
for nE=1:length(labels)
    for nG=1:length(uniqueGroupIDs)
        temp_plot=squeeze(all_distrib_SW_per(nE,:,all_GroupID==uniqueGroupIDs(nG)));
        simpleTplot(bins,temp_plot',0,Colors(nG,:),[0],'-',0.5,1,0,1,4);
    end
end
format_fig;
xlabel('Amplitude all waves')
ylabel('SW density')
xlim([0 75])

figure('Position',[347   308   441   306]);
uniqueGroupIDs={'Control','ADHD'};
hplot=[];
% for nG=1:length(uniqueGroupIDs)
for nE=1:length(labels)
    temp_plotA=squeeze(all_distrib_SW_per(nE,:,all_GroupID==uniqueGroupIDs(2)));
     temp_plotB=squeeze(all_distrib_SW_per(nE,:,all_GroupID==uniqueGroupIDs(1)));
     temp_plotC=squeeze(all_distrib_SW_per(nE,:,:));
     temp_plotdiff=temp_plotA-repmat(nanmean(temp_plotB,2),1,size(temp_plotA,2));

     distributionAll=cumsum(mean(temp_plotC,1))./max(cumsum(mean(temp_plotC,1)));

     bins_over=bins(find(mean(temp_plotA')>mean(temp_plotB')));
     bins_over(bins_over==0)=[];
     first_bin_over=bins_over(1);
     all_first_bin_over(nE)=first_bin_over;
%      all_first_bin_over(nE)=find(first_bin_over==bins)/length(bins);
     all_first_pctile_over(nE)=distributionAll(first_bin_over);
%    simpleTplot(bins,temp_plotA',0,Colors(2,:),[0 0.05 0.05 1000],'-',0.5,1,0,1,4);
hold on;
plot(bins,mean(temp_plotdiff'))
end
% legend(hplot,uniqueGroupIDs)
format_fig;
xlabel('Amplitude all waves')
ylabel('Diff. SW density')
xlim([0 75])
ylim([-2 2]*3)

%%
figure;
uniqueMSs={'ON','MW','MB','DK'};
cmap=cbrewer('qual','Set2',6);
hplot=[];
for nG=1:length(uniqueMSs)
    temp_plot=all_distrib_SW(SW_table.Probe_MS==uniqueMSs(nG),:);
    [pV hplot(nG)]=simpleTplot(bins,temp_plot,0,cmap(nG,:),[0],'-',0.5,1,0,1,4);
end
legend(hplot,uniqueMSs)
format_fig;
xlabel('Amplitude all waves')
ylabel('SW density')
xlim([0 100])

%%
figure;
cmap=cbrewer('seq','Blues',6);
hplot=[];
for nV=1:4
    temp_plot=all_distrib_SW(SW_table.Probe_Vig==nV,:);
    [pV hplot(nG)]=simpleTplot(bins,temp_plot,0,cmap(2+nV,:),[0],'-',0.5,1,0,1,4);
end
format_fig;
xlabel('Amplitude all waves')
ylabel('SW density')
xlim([0 100])

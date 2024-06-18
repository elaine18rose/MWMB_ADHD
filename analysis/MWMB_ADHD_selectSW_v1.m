%% Init
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
RS = ["R1", "R2"];
redo=0;
all_threshold_SW=readtable([preproc_path filesep 'all_threshold_SW.csv']);
all_threshold_SW.SubID=categorical(all_threshold_SW.SubID);
all_threshold_SW.Group=categorical(all_threshold_SW.Group);
all_threshold_SW.Elec=categorical(all_threshold_SW.Elec);
CTR_threshold_SW=all_threshold_SW(all_threshold_SW.SubID~='C017',:);

[Elec_group, ~, idx] = unique(CTR_threshold_SW.Elec);
mean_Thr = splitapply(@mean, CTR_threshold_SW.Thr_EG, idx);
av_CTR_threshold_SW = table(Elec_group, mean_Thr, 'VariableNames', {'Elec', 'mean_Thr'});

cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.channel=cellstr(Elec_group);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

figure;
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(av_CTR_threshold_SW.mean_Thr(av_CTR_threshold_SW.Elec==layout.label{nCh})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;

%%
SW_table=array2table(zeros(0,17),'VariableNames',{'SubID','Group','Block','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold','SW_peakneg','SW_peakpos','Probe_MS','Probe_Vig','Behav_Miss','Behav_FA','Behav_RT'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Group=categorical(SW_table.Group);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.Probe_MS=categorical(SW_table.Probe_MS);

MS_labels={'ON','MW','MB','DK'};
FilesPbme=[];
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
        FilesPbme=[FilesPbme ; {SubID} , {'Noisy SW detection'}];
        continue;
    end
    if strcmp(SubID,'C015')
        FilesPbme=[FilesPbme ; {SubID} , {'Weird behaviour'}];
        continue;
    end

    behav_file=dir([data_path filesep '..' filesep 'Behaviour' filesep 'wanderIM_behavres_' SubID '*.mat']);
    if length(behav_file)
        load([behav_file.folder filesep behav_file.name])
    else
        FilesPbme=[FilesPbme ; {SubID} , {'Missing Behaviour'}];
        continue;
    end

    fsample=500;
    labels={'Fp1'	'Fp2'	'F7'	'F3'	'Fz'	'F4'	'F8'	'FC5'	'FC1'	'FC2'	'FC6'	'T7'	'C3'	'Cz'	'C4'	'T8'	'TP9'	'CP5'	'CP1'	'CP2'	'CP6'	'TP10'	'P7'	'P3'	'Pz'	'P4'	'P8'	'PO9'	'O1'	'Oz'	'O2'	'PO10'	'AF7'	'AF3'	'AF4'	'AF8'	'F5'	'F1'	'F2'	'F6'	'FT9'	'FT7'	'FC3'	'FC4'	'FT8'	'FT10'	'C5'	'C1'	'C2'	'C6'	'TP7'	'CP3'	'CPz'	'CP4'	'TP8'	'P5'	'P1'	'P2'	'P6'	'PO7'	'PO3'	'POz'	'PO4'	'PO8'};
    load([preproc_path filesep 'SW_clean_i_probe_' SubID]);

    if size(probe_res,1)~=length(unique(all_Waves(:,2)))
        FilesPbme=[FilesPbme ; {SubID} , {'Different Numbers of Blocks'}];
        continue;
    end

    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
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

    thr_Wave_pc=[];
    thr_Wave_eg=[];
    slow_Waves=[];
    for nE=1:length(labels)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);

        this_thr=av_CTR_threshold_SW.mean_Thr(av_CTR_threshold_SW.Elec==labels{nE});

%         [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_p2p); % Fits an ex-Gauss distribution to data
%         bins=0:0.1:paramSW.art_ampl;                            % Creating variable to be used for the function below; from 0 to value set above (paramSW.art_ampl) in increments of 0.1
%         eg_pdf=exgauss_pdf(bins,X);                             % Computes the probability density; "X" here gives the values for [Mu, Sigma, Tau]
%         end_gaussian=2*bins(find(eg_pdf==max(eg_pdf)));
%         this_thr_ind=end_gaussian;
%         thr_Wave_eg(nE)=this_thr_ind;
        slow_Waves=[slow_Waves ; thisE_Waves(thisE_Waves(:,paramSW.AmpCriterionIdx)>this_thr,:)];
    end
%     save([preproc_path filesep 'selectSW_clean_i_probe_' SubID],'slow_Waves')

    for nBl=unique(slow_Waves(:,2))'
        slow_Waves_perE=[];
        duration_of_recording=20/60;
        for nE=1:length(labels)
            this_thr=av_CTR_threshold_SW.mean_Thr(av_CTR_threshold_SW.Elec==labels{nE});
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl)/duration_of_recording nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,5))/fsample)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,13)) this_thr nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,9)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,11))]];
        end

        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(labels)))=repmat({SubID},length(labels),1);
        SW_table.Group(table_length+(1:length(labels)))=repmat({GroupID},length(labels),1);
        SW_table.Block(table_length+(1:length(labels)))=repmat(nBl,length(labels),1);
        SW_table.Elec(table_length+(1:length(labels)))=labels;
        SW_table.SW_density(table_length+(1:length(labels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(labels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(labels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(labels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(labels)))=slow_Waves_perE(:,5);
        SW_table.SW_threshold(table_length+(1:length(labels)))=slow_Waves_perE(:,6);
        SW_table.SW_peakneg(table_length+(1:length(labels)))=slow_Waves_perE(:,7);
        SW_table.SW_peakpos(table_length+(1:length(labels)))=slow_Waves_perE(:,8);

        this_MS=MS_labels(probe_res(nBl,19));
        SW_table.Probe_MS(table_length+(1:length(labels)))=repmat(this_MS,length(labels),1);
        SW_table.Probe_Vig(table_length+(1:length(labels)))=repmat(probe_res(nBl,22),length(labels),1);

        this_block=probe_res(nBl,4);
        this_trial=probe_res(nBl,6);
        temp_gos=test_res(test_res(:,1)==this_block & test_res(:,4)<this_trial & test_res(:,5)~=3,:);
        temp_gos=temp_gos(end-17:end,:);
        temp_nogos=test_res(test_res(:,1)==this_block & test_res(:,4)<this_trial & test_res(:,5)==3,:);
        temp_nogos=temp_nogos(end-1:end,:);
        SW_table.Behav_Miss(table_length+(1:length(labels)))=repmat(nanmean(temp_gos(:,12)==0),length(labels),1);
        SW_table.Behav_RT(table_length+(1:length(labels)))=repmat(nanmean(temp_gos(:,10)-temp_gos(:,8)),length(labels),1);
        SW_table.Behav_FA(table_length+(1:length(labels)))=repmat(nanmean(temp_nogos(:,11)==0),length(labels),1);
    end

end
writetable(SW_table,[preproc_path filesep 'all_SW_perProbe_exGaussCTR.csv'])

%%
figure;
subplot(2,1,1)
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([0 16])

subplot(2,1,2)
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([0 16])
%%
figure;
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_threshold(SW_table.Elec==layout.label{nCh})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;

%%
figure;
for nState=1:4
    subplot(2,4,nState)
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Probe_MS==MS_labels{nState})));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    colormap(cmap); colorbar;
    caxis([0 16])
    title(MS_labels{nState})
end

for nVig=1:4
    subplot(2,4,4+nVig)
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Probe_Vig==nVig)));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    colormap(cmap); colorbar;
    caxis([0 16])
    title(sprintf('Vig%g',nVig))
end
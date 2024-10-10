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
RS = ["R1", "R2"];
%redo=0;
all_threshold_SW=readtable([preproc_path filesep 'all_threshold_SW_v2.csv']); %EP - This is
all_threshold_SW.SubID=categorical(all_threshold_SW.SubID);
all_threshold_SW.Group=categorical(all_threshold_SW.Group);
all_threshold_SW.Elec=categorical(all_threshold_SW.Elec);
CTR_threshold_SW=all_threshold_SW(all_threshold_SW.Group=='Control' & all_threshold_SW.SubID~='C017',:);

CTR_threshold_SW.Thr_EG(CTR_threshold_SW.Thr_EG>(nanmean(CTR_threshold_SW.Thr_EG)+5*nanstd(CTR_threshold_SW.Thr_EG)))=NaN;
[Elec_group, ~, idx] = unique(CTR_threshold_SW.Elec);
mean_Thr = splitapply(@nanmean, CTR_threshold_SW.Thr_EG, idx);
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
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colormap(cmap); colorbar;


load([pwd filesep '..' filesep 'Preproc' filesep 'all_badChannels_badProbes.mat']);

%%
SW_table=array2table(zeros(0,17),'VariableNames',{'SubID','Group','Block','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold','SW_peakneg','SW_peakpos','Probe_MS','Probe_Vig','Behav_Miss','Behav_FA','Behav_RT'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Group=categorical(SW_table.Group);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.Probe_MS=categorical(SW_table.Probe_MS);

SW_table=array2table(zeros(0,19),'VariableNames',{'SubID','Group','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold','SW_peakneg','SW_peakpos','Rate_ON','Rate_MW','Rate_MB','Rate_DR','Av_Vig','Behav_Miss','Behav_FA','Behav_RT'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Group=categorical(SW_table.Group);
SW_table.Elec=categorical(SW_table.Elec);

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

    if strcmp(SubID,'C017') || strcmp(SubID,'C038')
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
    %     labels={'Fp1'	'Fp2'	'F7'	'F3'	'Fz'	'F4'	'F8'	'FC5'	'FC1'	'FC2'	'FC6'	'T7'	'C3'	'Cz'	'C4'	'T8'	'TP9'	'CP5'	'CP1'	'CP2'	'CP6'	'TP10'	'P7'	'P3'	'Pz'	'P4'	'P8'	'PO9'	'O1'	'Oz'	'O2'	'PO10'	'AF7'	'AF3'	'AF4'	'AF8'	'F5'	'F1'	'F2'	'F6'	'FT9'	'FT7'	'FC3'	'FC4'	'FT8'	'FT10'	'C5'	'C1'	'C2'	'C6'	'TP7'	'CP3'	'CPz'	'CP4'	'TP8'	'P5'	'P1'	'P2'	'P6'	'PO7'	'PO3'	'POz'	'PO4'	'PO8'};
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
    if size(probe_res,1)~=length(unique(all_Waves(:,2)))
        FilesPbme=[FilesPbme ; {SubID} , {'Different Numbers of Probes'}];
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
        %
        %                 [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_p2p); % Fits an ex-Gauss distribution to data
        %                 bins=0:0.1:paramSW.art_ampl;                            % Creating variable to be used for the function below; from 0 to value set above (paramSW.art_ampl) in increments of 0.1
        %                 eg_pdf=exgauss_pdf(bins,X);                             % Computes the probability density; "X" here gives the values for [Mu, Sigma, Tau]
        %                 end_gaussian=2*bins(find(eg_pdf==max(eg_pdf)));
        %                 this_thr=end_gaussian;
        thr_Wave_eg(nE)=this_thr;
        slow_Waves=[slow_Waves ; thisE_Waves(thisE_Waves(:,paramSW.AmpCriterionIdx)>this_thr,:)];
    end
    %     save([preproc_path filesep 'selectSW_clean_i_probe_' SubID],'slow_Waves')

    for nBl=unique(slow_Waves(:,2))'
        slow_Waves_perE=[];
        duration_of_recording=20/60;
        for nE=1:length(labels)
            this_thr=thr_Wave_eg(nE); %av_CTR_threshold_SW.mean_Thr(av_CTR_threshold_SW.Elec==labels{nE});
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
        SW_table.Probe_ON(table_length+(1:length(labels)))=repmat(probe_res(nBl,19)==1,length(labels),1);
        SW_table.Probe_MW(table_length+(1:length(labels)))=repmat(probe_res(nBl,19)==2,length(labels),1);
        SW_table.Probe_MB(table_length+(1:length(labels)))=repmat(probe_res(nBl,19)==3,length(labels),1);
        SW_table.Probe_DR(table_length+(1:length(labels)))=repmat(probe_res(nBl,19)==4,length(labels),1);

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

    slow_Waves_perE=[];
    duration_of_recording=20/60*40;
    for nE=1:length(labels)
        this_thr=thr_Wave_eg(nE); %av_CTR_threshold_SW.mean_Thr(av_CTR_threshold_SW.Elec==labels{nE});
        slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE)/duration_of_recording nanmean(slow_Waves(slow_Waves(:,3)==nE,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE,7)-slow_Waves(slow_Waves(:,3)==nE,5))/fsample)) ...
            nanmean(slow_Waves(slow_Waves(:,3)==nE,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE,13)) this_thr nanmean(slow_Waves(slow_Waves(:,3)==nE,9)) nanmean(slow_Waves(slow_Waves(:,3)==nE,11))]];
    end

    table_length=size(SW_table_perS,1);
    SW_table_perS.SubID(table_length+(1:length(labels)))=repmat({SubID},length(labels),1);
    SW_table_perS.Group(table_length+(1:length(labels)))=repmat({GroupID},length(labels),1);
    SW_table_perS.Elec(table_length+(1:length(labels)))=labels;
    SW_table_perS.SW_density(table_length+(1:length(labels)))=slow_Waves_perE(:,1);
    SW_table_perS.SW_amplitude(table_length+(1:length(labels)))=slow_Waves_perE(:,2);
    SW_table_perS.SW_frequency(table_length+(1:length(labels)))=slow_Waves_perE(:,3);
    SW_table_perS.SW_downslope(table_length+(1:length(labels)))=slow_Waves_perE(:,4);
    SW_table_perS.SW_upslope(table_length+(1:length(labels)))=slow_Waves_perE(:,5);
    SW_table_perS.SW_threshold(table_length+(1:length(labels)))=slow_Waves_perE(:,6);
    SW_table_perS.SW_peakneg(table_length+(1:length(labels)))=slow_Waves_perE(:,7);
    SW_table_perS.SW_peakpos(table_length+(1:length(labels)))=slow_Waves_perE(:,8);

    SW_table_perS.Rate_ON(table_length+(1:length(labels)))=repmat(nanmean(probe_res(:,19)==1),length(labels),1);
    SW_table_perS.Rate_MW(table_length+(1:length(labels)))=repmat(nanmean(probe_res(:,19)==2),length(labels),1);
    SW_table_perS.Rate_MB(table_length+(1:length(labels)))=repmat(nanmean(probe_res(:,19)==3),length(labels),1);
    SW_table_perS.Rate_DR(table_length+(1:length(labels)))=repmat(nanmean(probe_res(:,19)==4),length(labels),1);
    SW_table_perS.Probe_Vig(table_length+(1:length(labels)))=repmat(nanmean(probe_res(:,22)),length(labels),1);


    temp_gos=test_res(test_res(:,5)~=3,:);
    temp_nogos=test_res(test_res(:,5)==3,:);
    SW_table_perS.Behav_Miss(table_length+(1:length(labels)))=repmat(nanmean(temp_gos(:,12)==0),length(labels),1);
    SW_table_perS.Behav_RT(table_length+(1:length(labels)))=repmat(nanmean(temp_gos(:,10)-temp_gos(:,8)),length(labels),1);
    SW_table_perS.Behav_FA(table_length+(1:length(labels)))=repmat(nanmean(temp_nogos(:,11)==0),length(labels),1);

end
% writetable(SW_table,[preproc_path filesep 'all_SW_perProbe_exGaussCTR_v2.csv'])

%%
figure;
subplot(2,1,1)
temp_topo=[];
temp_topo2=[];
for nCh=1:length(layout.label)-2
    temp_SW=SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD');
    temp_Sub=SW_table.SubID(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD');
    temp_av=grpstats(temp_SW(~isnan(temp_SW)),temp_Sub(~isnan(temp_SW)));
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD')));
    temp_topo2(nCh,:)=squeeze((SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([0 16])
title('ADHD')

subplot(2,1,2)
temp_topo=[];
temp_topo3=[];
for nCh=1:length(layout.label)-2
    temp_SW=SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control');
    temp_Sub=SW_table.SubID(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control');
    temp_av=grpstats(temp_SW(~isnan(temp_SW)),temp_Sub(~isnan(temp_SW)));
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control')));
    temp_topo3(nCh,:)=squeeze((SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([0 16])
title('Control')
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

figure;
subplot(2,1,1)
temp_topo=[];
temp_topo2=[];
for nCh=1:length(layout.label)-2
    temp_SW=SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD');
    temp_Sub=SW_table.SubID(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD');
    temp_av=grpstats(temp_SW(~isnan(temp_SW)),temp_Sub(~isnan(temp_SW)));
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD')));
    temp_topo2(nCh,:)=squeeze((SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='ADHD')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
% caxis([0 16])
title('Amplitude avADHD')

subplot(2,1,2)
temp_topo=[];
temp_topo3=[];
for nCh=1:length(layout.label)-2
    temp_SW=SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control');
    temp_Sub=SW_table.SubID(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control');
    temp_av=grpstats(temp_SW(~isnan(temp_SW)),temp_Sub(~isnan(temp_SW)));
    temp_topo(nCh)=squeeze(nanmean(SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control')));
    temp_topo3(nCh,:)=squeeze((SW_table.SW_amplitude(SW_table.Elec==layout.label{nCh} & SW_table.Group=='Control')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
% caxis([0 16])
title('Amplitude Control')

%% Ctr vs ADHD SW figure
f2=figure;
clear temp*

Groups={'ADHD','Control'};
VOI={'SW_density'}; %,'SW_amplitude','SW_frequency','SW_downslope','SW_upslope'};
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

for nP=1:length(VOI)
    minmax_val=[];
    for nGroup=1:2
        figure(f2);
        sgtitle('Control v ADHD')
        subplot(2,length(VOI),nP+(nGroup-1)*length(VOI))


        temp_topo=[];
        for nE=1:length(layout.label)-2
            temp=SW_table.(VOI{nP})(SW_table.Elec==layout.label{nE} & SW_table.Group==Groups{nGroup}); %Getting the values per electrode by group
            temp_topo(nE)=nanmean(temp); % Getting the means per electrode by group
        end

        simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
        colorbar;
        cmap=colormap('hot'); cmap=flipud(cmap);
        colormap(cmap);
        title({VOI{nP}(4:end),Groups{nGroup}})
        minmax_val=[minmax_val ; min(temp_topo) max(temp_topo)];
        caxis([5 18])
        format_fig;
    end

end

%% Topos of the relationship between Group and SWDens using LMEs with BlockN as random slope
clear topo_*
clear sub_table
clear temp_mdl_*
SW_table.Group=categorical(SW_table.Group);
SW_table.Group=reordercats(SW_table.Group,["Control" "ADHD"]);

for nE=1:length(layout.label)-2
    fprintf('... %g/%g\n',nE,length(layout.label)-2)
    sub_table=SW_table(SW_table.Elec==layout.label(nE),:);
    temp_mdl_byG=fitlme(sub_table,sprintf('SW_density~1+%s+Block+(1|SubID)','Group'));

    topo_tV_byGroup(nE)=temp_mdl_byG.Coefficients.tStat(match_str(temp_mdl_byG.Coefficients.Name,"Group_ADHD"));
    topo_pV_byGroup(nE)=temp_mdl_byG.Coefficients.pValue(match_str(temp_mdl_byG.Coefficients.Name,"Group_ADHD"));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
figure;
simpleTopoPlot_ft(topo_tV_byGroup', layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(topo_pV_byGroup<0.05), 'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',72,'box','no','label','no');
%ft_plot_lay_me(layout, 'chanindx', find(topo_pV_byGroup<fdr(topo_pV_byGroup,0.05)), 'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',72,'box','no','label','no');
colormap(cmap2);
colorbar('Ticks',[-5:5]);
title({'SW density','A vs C'})
caxis([-1 1]*4)
format_fig;

%%
SWdens_est=cell(1,2);
totperm=500;
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label(nCh),:);
    if nCh==1
        out_pred_perm=[];
        [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_bygroup(sub_table,'Group','SW_density~1+Block+pred+(1|SubID)',totperm);
    else
        [real_out, cont_out, perm_out, cont_perm_out, next_out_pred_perm]=lme_perm_bygroup(sub_table,'Group','SW_density~1+Block+pred+(1|SubID)',totperm,out_pred_perm);
    end
    SWdens_est{1}=[SWdens_est{1} ; [nCh real_out]];
    SWdens_est{2}=[SWdens_est{2} ; [nCh*ones(totperm,1) perm_out]];
end

%%
clus_alpha=0.1;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='EEG1010.lay';
cfg_neighb.channel=layout.label(1:end-2);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(SW_table.Elec)))=[];
[SWdens_clus]=get_clusterperm_lme_bygroup(SWdens_est,clus_alpha,montecarlo_alpha,totperm,neighbours,1);


cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=1;
limMax=5;

temp_topo=SWdens_est{1}(:,3);
temp_topo2=zeros(size(temp_topo));
temp_topo3=zeros(size(temp_topo));
temp_clus=SWdens_clus;

figure;
if ~isempty(temp_clus)
    for nclus=1:length(temp_clus)
        if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
            continue;
        end
        %             ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','yes')
        temp_topo2(match_str(layout.label,temp_clus{nclus}{2}))=temp_topo(match_str(layout.label,temp_clus{nclus}{2}));
        temp_topo3(match_str(layout.label,temp_clus{nclus}{2}))=1;
        fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
    end
end
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
%     ft_plot_lay_me(layout, 'chanindx',1:length(layout.label)-2,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
format_fig;
caxis([-1 1]*limMax)
if ~isempty(temp_clus)
    for nclus=1:length(temp_clus)
        if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
            continue;
        end
        ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',96,'box','no','label','no')
    end
end
hb=colorbar('Position',[0.94    0.6    0.05    0.33]);
colormap(cmap2);

%% Correlation with behaviour
clear topo_*
clear sub_table
clear temp_mdl_*
SW_table.Group=categorical(SW_table.Group);
SW_table.Group=reordercats(SW_table.Group,["Control" "ADHD"]);
VOI={'Behav_FA','Behav_Miss','Behav_RT','Probe_Vig','Probe_MW','Probe_MB'};
figure;
for nV=1:length(VOI)
    topo_tV=[];
    topo_pV=[];
    for nE=1:length(layout.label)-2
        fprintf('... %g/%g\n',nE,length(layout.label)-2)
        sub_table=SW_table(SW_table.Elec==layout.label(nE),:);
        temp_mdl_byBehav=fitlme(sub_table,sprintf('SW_density~1+Group+%s+Block+(1|SubID)',VOI{nV}));

        topo_tV(nE)=temp_mdl_byBehav.Coefficients.tStat(find_trials(temp_mdl_byBehav.Coefficients.Name,VOI{nV}));
        topo_pV(nE)=temp_mdl_byBehav.Coefficients.pValue(find_trials(temp_mdl_byBehav.Coefficients.Name,VOI{nV}));
    end
    subplot(2,3,nV)
    cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
    simpleTopoPlot_ft(topo_tV', layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx', find(topo_pV<0.05), 'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',72,'box','no','label','no');
    %ft_plot_lay_me(layout, 'chanindx', find(topo_pV_byGroup<fdr(topo_pV_byGroup,0.05)), 'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',72,'box','no','label','no');
    colormap(cmap2);
    colorbar('Ticks',[-5:5]);
    title(VOI{nV})
    caxis([-1 1]*max(abs(topo_tV)))
    format_fig;
end
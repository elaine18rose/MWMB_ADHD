%% Copied from CTET ADHD Script 
clear all
close all

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

cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

%%
nFc=0;
thr_Wave=[];
 all_slow_Waves=[];
    all_SubIDs=[];
all_Groups=[];
for nF=1:length(SW_files)
    file_name = SW_files(nF).name;
    folder_name = SW_files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(seps(end)+1:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(SW_files))
    
    behav_file=dir([data_path filesep '..' filesep 'Behaviour' filesep 'wanderIM_behavres_' SubID '*.mat']);
    if length(behav_file)
        load([behav_file.folder filesep behav_file.name])
    else
        FilesPbme=[FilesPbme ; {SubID} , {'Missing Behaviour'}];
        continue;
    end
    
    load([preproc_path filesep 'SW_clean_i_probe_' SubID]); %,'all_Waves')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 8]; % [1 4] delta or [4 8] theta
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
%     paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    fsample=256;
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2))*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    slow_Waves=[];
    nFc=nFc+1;
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
%         if ~isempty(paramSW.fixThr)
%             thr_Wave(nFc,nE)=paramSW.fixThr;
%         else
%             thr_Wave(nFc,nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
%         end
    
        slow_Waves=[slow_Waves ; thisE_Waves];
    end
    save([preproc_path filesep filesep 'noThr_SW_clean_i_probe_' SubID],'slow_Waves','paramSW')
    
    all_slow_Waves=[ all_slow_Waves ; slow_Waves ];
    all_SubIDs=[all_SubIDs ; repmat({SubID},size(slow_Waves,1),1)];
    
    orifoldername=SW_files(nF).name;
    if isempty(findstr(orifoldername,'SW_clean_i_probe_C'))==0
        group_SW={'Control'};
    elseif isempty(findstr(orifoldername,'SW_clean_i_probe_A'))==0
        group_SW={'ADHD'};
    end
    all_Groups=[all_Groups ; repmat(group_SW,size(slow_Waves,1),1)];
end

%%
% Col  1: subject
                    % Col  2: epoch number
                    % Col  3: electrode number (as in data_clean.labels) (***)
                    % Col  4: peak to peak amplitude (***)
                    % Col  5: start wave (sample same sampling rate as data_clean) (***)
                    % Col  6: mid wave
                    % Col  7: end wave
                    % Col  8: negative peak position
                    % Col  9: negative peak amplitude
                    % Col 10: positive peak position
                    % Col 11: positive peak amplitude
                    % Col 12: downward slope
                    % Col 13: upward slope
                    % Col 14: max amplitude around wave
                    % Col 15: min amplitude around wave
                    
all_SW_table=array2table(all_slow_Waves,'VariableNames',{'SubID','ProbeN','ElecN','Amplitude','Start','Middle','End','NegPeak','NegPeakAmp','PosPeak','PosPeakAmp','Down_Slope','Up_Slope','MaxAmp','MinAmp'});
all_SW_table.SubID=all_SubIDs;
all_SW_table.Group=all_Groups;

writetable(all_SW_table,[preproc_path filesep 'noThr_SW_clean_i_probe_All.csv']);

%%
SubIDs=unique(all_SW_table.SubID);
Probes=unique(all_SW_table.ProbeN);
Elecs=unique(all_SW_table.ElecN);

av_SW_matrix=nan(length(SubIDs)*length(Probes)*length(Elecs),9);
all_SubID=cell(length(SubIDs)*length(Probes)*length(Elecs),1);
all_Group=cell(length(SubIDs)*length(Probes)*length(Elecs),1);
count=0;
for nF=1:length(SubIDs)
    for nProbe=1:length(Probes)
                    beg_block=min(all_SW_table.Start(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ProbeN,Probes(nProbe))));
            end_block=max(all_SW_table.End(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ProbeN,Probes(nProbe))));
            duration_block=(end_block-beg_block)/fsample/60;
            
        for nElec=1:length(Elecs)
            temp_table=all_SW_table(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ProbeN,Probes(nProbe)) & ismember(all_SW_table.ElecN,Elecs(nElec)),:);
         
            SW_dens=size(temp_table,1);

            SW_dens=SW_dens/duration_block;

            count=count+1;
            if isempty(SW_dens)
            av_SW_matrix(count,:)=[nF nProbe nElec 0 nan(1,5)];
            else
            av_SW_matrix(count,:)=[nF nProbe nElec SW_dens mean(table2array(temp_table(:,[4 9 11 12 13])),1)];
            end
            temp_table2=all_SW_table(ismember(all_SW_table.SubID,SubIDs(nF)),:);
            all_SubID(count)=SubIDs(nF);
            all_Group(count)=unique(temp_table2.Group);
            fprintf('%g - %g - %g\n',nF,nProbe,nElec)
            if size(av_SW_matrix,1)~=length(all_SubID)
                warning('oops')
                pause;
            end
        end
    end
end
av_SW_table=array2table(av_SW_matrix,'VariableNames',{'SubID','nProbe','ElecN','SWdens','Amplitude','NegPeakAmp','PosPeakAmp','Down_Slope','Up_Slope'});

av_SW_table.Group=all_Group;
av_SW_table.SubID=all_SubID;

writetable(av_SW_table,[preproc_path filesep  'noThr_SW_clean_i_probe_perProbe.csv']);

%%
SubIDs=unique(all_SW_table.SubID);
Probes=unique(all_SW_table.ProbeN);
Elecs=unique(all_SW_table.ElecN);

av_SW_matrix=nan(length(SubIDs)*length(Elecs),9);
all_SubID=cell(length(SubIDs)*length(Elecs),1);
all_Group=cell(length(SubIDs)*length(Elecs),1);
count=0;
for nF=1:length(SubIDs)
    duration_probe=0;
    for nProbe=1:length(Probes)
        beg_block=min(all_SW_table.Start(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ProbeN,Probes(nProbe))));
        end_block=max(all_SW_table.End(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ProbeN,Probes(nProbe))));
        duration_probe=duration_probe + (end_block-beg_block)/fsample/60;
    end

        for nElec=1:length(Elecs)
            tic;
            temp_table=all_SW_table(ismember(all_SW_table.SubID,SubIDs(nF)) & ismember(all_SW_table.ElecN,Elecs(nElec)),:);
         
            SW_dens=size(temp_table,1);

            SW_dens=SW_dens/duration_probe;

            count=count+1;
            if isempty(SW_dens)
            av_SW_matrix(count,:)=[nF 0 nElec 0 nan(1,5)];
            else
            av_SW_matrix(count,:)=[nF 0 nElec SW_dens mean(table2array(temp_table(:,[4 9 11 12 13])),1)];
            end
            temp_table2=all_SW_table(ismember(all_SW_table.SubID,SubIDs(nF)),:);
            all_SubID(count)=SubIDs(nF);
            all_Group(count)=unique(temp_table2.Group);
            toc;
            fprintf('%g 0- %g\n',nF,nElec)
            if size(av_SW_matrix,1)~=length(all_SubID)
                warning('oops')
                pause;
            end
        end
end
av_SW_table=array2table(av_SW_matrix,'VariableNames',{'SubID','ProbeN','ElecN','SWdens','Amplitude','NegPeakAmp','PosPeakAmp','Down_Slope','Up_Slope'});

av_SW_table.Group=all_Group;
av_SW_table.SubID=all_SubID;

writetable(av_SW_table,[preproc_path filesep  'noThr_SW_clean_i_probe.csv']);


%% Relabeling Elec No. to Labels in the table with average SWs
% Convert the 'ElecN' column to a cell array
av_SW_table.ElecN = num2cell(av_SW_table.ElecN);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for i = 1:height(av_SW_table)
    elec = av_SW_table.ElecN{i}; % Current electrode number (now a cell)
    if elec > 0 && elec <= length(labels)
        av_SW_table.ElecN{i} = labels{elec}; % Replace the number with the corresponding label
    else
        warning('Electrode number %d is out of range', elec);
    end
end

%% Figures 
% Swdensity
figure;
subplot(2,1,1)
temp_topo=[];
for nCh=1:length(labels)'
    temp_topo(nCh)=squeeze(nanmean(av_SW_table.SWdens(strcmp(labels{nCh}, av_SW_table.ElecN) & strcmp(av_SW_table.Group, 'ADHD'))));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([30 45])
title('ADHD - SWD')

subplot(2,1,2)
temp_topo=[];
for nCh=1:length(labels)'
    temp_topo(nCh)=squeeze(nanmean(av_SW_table.SWdens(strcmp(labels{nCh}, av_SW_table.ElecN) & strcmp(av_SW_table.Group, 'Control'))));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([30 45])
title('Control - SWD')


% Amp
figure;
subplot(2,1,1)
temp_topo=[];
for nCh=1:length(labels)'
    temp_topo(nCh)=squeeze(nanmean(av_SW_table.Amplitude(strcmp(labels{nCh}, av_SW_table.ElecN) & strcmp(av_SW_table.Group, 'ADHD'))));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([8 20])
title('ADHD - Amp')

subplot(2,1,2)
temp_topo=[];
for nCh=1:length(labels)'
    temp_topo(nCh)=squeeze(nanmean(av_SW_table.Amplitude(strcmp(labels{nCh}, av_SW_table.ElecN) & strcmp(av_SW_table.Group, 'Control'))));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
caxis([8 20])
title('Control - Amp')


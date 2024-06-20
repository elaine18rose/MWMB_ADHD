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
all_badCompo=[];
redo=1;
all_threshold_SW=array2table(zeros(0,5),'VariableNames',{'SubID','Group','Elec','Thr_PC','Thr_EG'});
all_threshold_SW.SubID=categorical(all_threshold_SW.SubID);
all_threshold_SW.Group=categorical(all_threshold_SW.Group);
all_threshold_SW.Elec=categorical(all_threshold_SW.Elec);
for nF=1:length(SW_files)
    if startsWith(SW_files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end

    if contains(SW_files(nF).name,RS) %To skip resting state files
        continue
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
    if redo==1 || exist([preproc_path filesep 'SW_clean_i_probe_' SubID '.mat'])==0 % To skip already preprocessed files
        load([preproc_path filesep 'clean_i_probe_' SubID '.mat']);

        all_Waves=[];
        for nBl=1:size(data.trial,2)
            temp_data=data.trial{nBl}(:,:);
            temp_data=temp_data-repmat(mean(temp_data(match_str(data.label,{'TP9','TP10'}),:),1),size(temp_data,1),1);
            temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));

            [twa_results]=twalldetectnew_TA_v2(temp_data,data.fsample,0);
            for nE=1:length(data.label)
                if ~isfield(twa_results.channels(nE),'maxnegpkamp')
                    continue;
                end
                all_Waves=[all_Waves ; [repmat([nF nBl nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                    cell2mat(twa_results.channels(nE).negzx)' ...
                    cell2mat(twa_results.channels(nE).poszx)' ...
                    cell2mat(twa_results.channels(nE).wvend)' ...
                    cell2mat(twa_results.channels(nE).maxnegpk)' ...
                    cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                    cell2mat(twa_results.channels(nE).maxpospk)' ...
                    cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                    cell2mat(twa_results.channels(nE).mxdnslp)' ...
                    cell2mat(twa_results.channels(nE).mxupslp)' ...
                    cell2mat(twa_results.channels(nE).maxampwn)' ...
                    cell2mat(twa_results.channels(nE).minampwn)' ...
                    ]];
            end
        end
        fprintf('\n')
        fsample=data.fsample;
        labels=data.label;
        save([preproc_path filesep 'SW_clean_i_probe_' SubID],'all_Waves','fsample','labels')
    else
        fsample=500;
        labels={'Fp1'	'Fp2'	'F7'	'F3'	'Fz'	'F4'	'F8'	'FC5'	'FC1'	'FC2'	'FC6'	'T7'	'C3'	'Cz'	'C4'	'T8'	'TP9'	'CP5'	'CP1'	'CP2'	'CP6'	'TP10'	'P7'	'P3'	'Pz'	'P4'	'P8'	'PO9'	'O1'	'Oz'	'O2'	'PO10'	'AF7'	'AF3'	'AF4'	'AF8'	'F5'	'F1'	'F2'	'F6'	'FT9'	'FT7'	'FC3'	'FC4'	'FT8'	'FT10'	'C5'	'C1'	'C2'	'C6'	'TP7'	'CP3'	'CPz'	'CP4'	'TP8'	'P5'	'P1'	'P2'	'P6'	'PO7'	'PO3'	'POz'	'PO4'	'PO8'};
        load([preproc_path filesep 'SW_clean_i_probe_' SubID]);
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

            thr_Wave_pc(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);

            [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_p2p); % Fits an ex-Gauss distribution to data
            bins=0:0.1:paramSW.art_ampl;                            % Creating variable to be used for the function below; from 0 to value set above (paramSW.art_ampl) in increments of 0.1
            eg_pdf=exgauss_pdf(bins,X);                             % Computes the probability density; "X" here gives the values for [Mu, Sigma, Tau]
            end_gaussian=2*bins(find(eg_pdf==max(eg_pdf)));         % Finding the end of the Gauss
            thr_Wave_eg(nE)=end_gaussian;      % Otherwise, if paramSW.fixThr is empty, save the value obtained in end_gaussian (.fixThr should be empty)
        end

     
        this_thr=array2table(zeros(length(labels),5),'VariableNames',{'SubID','Group','Elec','Thr_PC','Thr_EG'});
        this_thr.SubID=categorical(this_thr.SubID);
        this_thr.Group=categorical(this_thr.Group);
        this_thr.Elec=categorical(this_thr.Elec);

        this_thr.SubID=repmat({SubID},length(labels),1);
        this_thr.Group=repmat({GroupID},length(labels),1);
        this_thr.Elec=labels';
        this_thr.Thr_PC=thr_Wave_pc';
        this_thr.Thr_EG=thr_Wave_eg';

        all_threshold_SW=[all_threshold_SW ; this_thr];
end

writetable(all_threshold_SW,[preproc_path filesep 'all_threshold_SW.csv'])


%%
figure;
simpleCorPlot(all_threshold_SW.Thr_PC,all_threshold_SW.Thr_EG);
%simpleCorPlot(all_threshold_SW.Thr_PC(all_threshold_SW.SubID~='C017'),all_threshold_SW.Thr_EG(all_threshold_SW.SubID~='C017'));
xlabel('Percentile threshold')
ylabel('Ex-Gaussian threshold')
%title('without C017')

%%
cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.channel=labels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);


figure;
subplot(2,1,1);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(all_threshold_SW.Thr_PC(all_threshold_SW.Elec==layout.label{nCh} & all_threshold_SW.SubID~='C017')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;

subplot(2,1,2);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(all_threshold_SW.Thr_EG(all_threshold_SW.Elec==layout.label{nCh} & all_threshold_SW.SubID~='C017')));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;


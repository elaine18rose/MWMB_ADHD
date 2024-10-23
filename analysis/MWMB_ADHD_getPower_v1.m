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
eeg_files=dir([data_path filesep '*.eeg']);

%EEG Layout info
run ../MWMB_ADHD_elec_layout.m


%% Loop across files
RS = ["R1", "R2"];
all_badCompo=[];
redo=0;
all_pow=[];
all_frac=[];
all_osci=[];
all_aperiodic=[];
all_alphapeak=[];
nFc=0;
for nF=1:length(eeg_files)
    if startsWith(eeg_files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end

    if contains(eeg_files(nF).name,RS) %To skip resting state files
        continue
    end
    nFc=nFc+1;
    %%% load the data
    SubInfo=split(eeg_files(nF).name,'-');
    SubID=SubInfo{2}(1:end-4);
    if SubID(1)=='A'
        GroupID='ADHD';
    elseif SubID(1)=='C'
        GroupID='Control';
    else
        GroupID='undefined';
    end
    if exist([preproc_path filesep 'clean_i_probe_' SubID '.mat'])==0
        continue;
    end
    if redo==1 || exist([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat'])==0 
        load([preproc_path filesep 'clean_i_probe_' SubID '.mat']);

        % take the 20s before a probe
        cfg=[];
        cfg.toilim    = [-20 0];
        data          = ft_redefinetrial(cfg, data);

        cfg                             = [];
        cfg.foilim                      = [1 40]; %[1 30];
        cfg.pad                         = 'nextpow2';
        cfg.tapsmofrq                   = 0.5;
        cfg.method                      = 'mtmfft';
        cfg.output                      = 'fooof_aperiodic';
        cfg.fooof.max_peaks             = 4;
        cfg.fooof.proximity_threshold   = 1;
        cfg.keeptrials                  = 'no';
        fractal = ft_freqanalysis(cfg, data);

        cfg.keeptrials                  = 'no';
        cfg.output                      = 'pow';
        pow = ft_freqanalysis(cfg, data);

        cfg.keeptrials                  = 'no';
        cfg.output                      = 'fooof_peaks';
        pow_peaks = ft_freqanalysis(cfg, data);

        save([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat'],'pow','fractal','pow_peaks');
    else
        load([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat']);
    end
    all_pow(nFc,:,:)=log10(pow.powspctrm);
    all_frac(nFc,:,:)=fractal.powspctrm;
    all_osci(nFc,:,:)=pow_peaks.powspctrm;


    aperiodic_params=[];
    for nE=1:length(fractal.label)
        aperiodic_params(nE,:)=fractal.fooofparams(nE).aperiodic_params;
    end
    all_aperiodic(nFc,:,:) = aperiodic_params;

    peaks_params=[];
    for nE=1:length(fractal.label)
        temp_peaks=fractal.fooofparams(nE).peak_params;
        temp_peaks=temp_peaks(temp_peaks(:,1)>9 & temp_peaks(:,1)<13,:);
        if size(temp_peaks,1)==1
            peaks_params(nE,:)=temp_peaks(temp_peaks(:,2)==max(temp_peaks(:,2)),:);
        else
            peaks_params(nE,:)=nan(1,3);
        end
    end
    all_alphapeak(nFc,:,:) = peaks_params;

    % Copied from CTET script: 
        if SubID(1)=='C'
        group_PowData{nFc}='Control';
        design_PowData(1,nFc)=0;
    elseif SubID(1)=='A'
        group_PowData{nFc}='ADHD';
        design_PowData(1,nFc)=1;
    end
end

%%
figure;
subplot(1,3,1);
plot(pow.freq,squeeze(mean(all_pow,1))');
format_fig;
xlabel('Freq (Hz)'); title('Full Power');

subplot(1,3,2);
plot(pow.freq,squeeze(mean(all_osci,1))');
format_fig;
xlabel('Freq (Hz)'); title('Periodic');

subplot(1,3,3);
plot(pow.freq,squeeze(mean(all_frac,1))');
format_fig;
xlabel('Freq (Hz)'); title('Aperiodic');

%%
cfg = [];
cfg.layout = 'EEG1010.lay';
%cfg.channel=data.label;
cfg.channel=pow.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

freq_bands=[1 4;4 8;8 12];
name_bands={'\delta','\theta','\alpha'};
for nF=1:size(freq_bands,1)
    temp_topo{nF}=[];
    for nCh=1:length(layout.label)-2
         temp_topo{nF}(nCh)=squeeze(nanmean(nanmean(all_pow(:,match_str(pow.label,layout.label{nCh}),pow.freq>freq_bands(nF,1) & pow.freq<freq_bands(nF,2)),3),1));
%         temp_topo{nF}(nCh)=squeeze(nanmean(nanmean(all_pow(:,match_str(data.label,layout.label{nCh}),pow.freq>freq_bands(nF,1) & pow.freq<freq_bands(nF,2)),3),1));
        %         temp_topo{nF}(nCh)=squeeze(nanmean(nanmean(all_osci(:,match_str(data.label,layout.label{nCh}),pow.freq>freq_bands(nF,1) & pow.freq<freq_bands(nF,2)),3),1));
    end
end

cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
figure;
for nF=1:size(freq_bands,1)
    subplot(1,size(freq_bands,1),nF);

    simpleTopoPlot_ft(temp_topo{nF}', layout,'on',[],0,1);
    colormap(cmap); colorbar;
    title(name_bands{nF});
    format_fig;
end

%%
figure;
load([preproc_path filesep 'clean_i_probe_' SubID '.mat']); % Uncomment if crashing in 192
subplot(1,2,1);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(all_aperiodic(:,match_str(data.label,layout.label{nCh}),1),1));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
title('offset');
format_fig;

subplot(1,2,2);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(nanmean(all_aperiodic(:,match_str(data.label,layout.label{nCh}),2),1));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colormap(cmap); colorbar;
title('slope');
format_fig;

%% PS by group 
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

faxis = pow.freq;
chLabels=pow.label;

figure;
subplot(2,2,1);
hp=[];
[~,hp(1)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'Control'),match_str(chLabels,'Cz'),:))),0,Colors(1,:),0,'-',0.2,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'ADHD'),match_str(chLabels,'Cz'),:))),0,Colors(2,:),0,'-',0.2,1,0,1,2);
hold on;
xlim([1 40]) %xlim([1 10])
legend(hp,{'Controls','ADHD'})
title('Power spectrum at electrode Cz');
xlabel('Frequency (Hz)')
ylabel('Log Pow')
format_fig;

subplot(2,2,2);
hp=[];
[~,hp(1)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'Control'),match_str(chLabels,'Fz'),:))),0,Colors(1,:),0,'-',0.2,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'ADHD'),match_str(chLabels,'Fz'),:))),0,Colors(2,:),0,'-',0.2,1,0,1,2);
hold on;
xlim([1 40])
legend(hp,{'Controls','ADHD'})
title('Power spectrum at electrode Fz');
xlabel('Frequency (Hz)')
ylabel('Log Pow')
format_fig;

subplot(2,2,3);
hp=[];
[~,hp(1)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'Control'),match_str(chLabels,'Pz'),:))),0,Colors(1,:),0,'-',0.2,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'ADHD'),match_str(chLabels,'Pz'),:))),0,Colors(2,:),0,'-',0.2,1,0,1,2);
hold on;
xlim([1 40])
legend(hp,{'Controls','ADHD'})
title('Power spectrum at electrode Pz');
xlabel('Frequency (Hz)')
ylabel('Log Pow')
format_fig;

subplot(2,2,4);
hp=[];
[~,hp(1)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'Control'),match_str(chLabels,'Oz'),:))),0,Colors(1,:),0,'-',0.2,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(faxis,squeeze((all_pow(match_str(group_PowData,'ADHD'),match_str(chLabels,'Oz'),:))),0,Colors(2,:),0,'-',0.2,1,0,1,2);
hold on;
xlim([1 40])
legend(hp,{'Controls','ADHD'})
title('Power spectrum at electrode Oz');
xlabel('Frequency (Hz)')
ylabel('Log Pow')
format_fig;


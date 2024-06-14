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
eeg_files=dir([data_path filesep '*.eeg']);

%EEG Layout info
run ../MWMB_ADHD_elec_layout.m


%% Loop across files
all_badCompo=[];
redo=1;
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
    elseif SubID(2)=='C'
        GroupID='Control';
    else
        GroupID='undefined';
    end
    if redo==1 || exist([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat'])==0 % To skip already preprocessed files
        load([preproc_path filesep 'clean_i_probe_' SubID '.mat']);

        % take the 20s before a probe
        cfg=[];
        cfg.toilim    = [-20 0];
        data          = ft_redefinetrial(cfg, data);
        
        cfg                             = [];
        cfg.foilim                      = [1 30];
        cfg.pad                         = 'nextpow2';
        cfg.tapsmofrq                   = 0.5;
        cfg.method                      = 'mtmfft';
        cfg.output                      = 'fooof_aperiodic';
        cfg.fooof.max_peaks             = 4;
        cfg.fooof.proximity_threshold   = 1;
        cfg.keeptrials                  = 'no';
        fractal = ft_freqanalysis(cfg, data);
        
        cfg.keeptrials                  = 'yes';
        cfg.output                      = 'pow';
        pow = ft_freqanalysis(cfg, data);
        
        cfg.keeptrials                  = 'no';
        cfg.output                      = 'fooof_peaks';
        pow_peaks = ft_freqanalysis(cfg, data);
                
        save([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat'],'pow','fractal','pow_peaks');
    else
        load([preproc_path filesep 'TF_clean_i_probe_' SubID '.mat']);
    end
    all_pow(nFc,:,:,:)=log10(pow.powspctrm);
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
    

end



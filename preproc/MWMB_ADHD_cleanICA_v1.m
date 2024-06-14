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
addpath(path_eeglab);
eeglab;

% select relevant files, here task
eeg_files=dir([data_path filesep '*.eeg']);

%EEG Layout info
run ../MWMB_ADHD_elec_layout.m


%% Loop across files
RS = ["R1", "R2"];
all_badCompo=[];
redo=1;
all_ICA_classification=[];
for nF=1:length(eeg_files)
    if startsWith(eeg_files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end

    if contains(eeg_files(nF).name,RS) %To skip resting state files
        continue
    end

    %%% load the data
    SubInfo=split(eeg_files(nF).name,'-');
    SubID=SubInfo{2}(1:end-4);

    if redo==1 || exist([preproc_path filesep 'clean_i_probe_' SubID '.mat'])==0 % To skip already preprocessed files
        if exist([preproc_path filesep 'comp_i_probe_' SubID '.mat'])==0
            warning(sprintf('... missing ICA file for %s\n',[eeg_files(nF).name]))
            continue;
        else
            fprintf('... working on %s\n',[eeg_files(nF).name])
        end

        load([preproc_path filesep 'comp_i_probe_' SubID])

        rejected_comps = ICA_classification.Comp(ICA_classification.Eye>0.95 | ICA_classification.Heart>0.8); % reject eye component with proba over 0.95. You could add heart compo as well.
        fprintf('... ... %g bad components rejected\n',length(rejected_comps))
        EEG_clean = pop_subcomp(EEG_ica, rejected_comps);
        EEG_clean = eeg_checkset(EEG_clean);
        badCompo=ICA_classification(rejected_comps,:);
        %             figure;
        %             plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'Fp1'),:,1))','r');
        %             hold on
        %             plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'Fp2'),:,1))','b');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Fp1'),:,1))','r--');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Fp2'),:,1))','b--');
        %
        %                 figure;
        %             plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'O1'),:,1))','r');
        %             hold on
        %             plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'O2'),:,1))','b');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'O1'),:,1))','r--');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'O2'),:,1))','b--');

        data = eeglab2fieldtrip( EEG_clean, 'raw');
        %  figure;
        %             plot(squeeze(data.trial{1}(match_str(data.label,'Fz'),:)),'r');
        %             hold on
        %             plot(squeeze(data.trial{1}(match_str(data.label,'Cz'),:)),'b');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Fz'),:,1)),'m--');
        %             plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Cz'),:,1)),'c--');
        %
        save([preproc_path filesep 'clean_i_probe_' SubID '.mat'],'data','badCompo');
        all_badCompo=[all_badCompo ; badCompo];
    end
end
writetable(all_badCompo,[preproc_path filesep 'ICA_classification_allSubs_badComponents.csv'])

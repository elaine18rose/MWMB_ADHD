%%
clear all;
close all;

%%
if isempty(findstr(pwd,'thandrillon'))==0
    preproc_path='/Users/thandrillon/Data/ADHD_MW/Preproc/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_data='/Users/thandrillon/Data/ADHD_MW/';
    path_eeglab='/Users/thandrillon/WorkGit/projects/ext/eeglab/';
    path_ICAlabel='/Users/thandrillon/WorkGit/projects/ext/ICLabel/';
else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/Desktop/MATLAB_Functions/fieldtrip/';
    path_data = '/Volumes/Seagate/MWMB_ADHD_SART/';

    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    path_detectSW = '/Volumes/Seagate/MWMB_ADHD_SART/SW_detection/';
    path_eeglab='/Users/elaine/desktop/MATLAB_Functions/eeglab/';
    path_ICAlabel='/Users/elaine/desktop/MATLAB_Functions/ICLabel/';

end

addpath(genpath(path_LSCPtools))
addpath(path_fieldtrip)
ft_defaults;
addpath(path_eeglab);
eeglab;

files=dir([path_data filesep 'EEG' filesep '*_R*.eeg']);

%% loop on subjects
redo=1; % 1 to force re-import, 0 otherwise
nFc=0;
all_ICA_classification=[];
all_badCompo=[];
for nF=1:length(files)
    if startsWith(files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    if isempty(findstr(SubID,'ID-'))==0
        SubID(1:findstr(SubID,'ID-')+2)=[];
    end
    FileNameID=SubID;
    if isempty(findstr(SubID,'_R'))==0
        SubID(findstr(SubID,'_R'):end)=[];
    end
    tic;
    %     if ~strcmp(SubID,'A084') %%re-add this if you need to clip the data
    %         continue;
    %     end
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))

    if redo==1 || exist([preproc_path filesep 'clean_pow_RS_' SubID '.mat'])==0

        %%% minimal preprocessing
        cfg=[];
        cfg.SubID          = SubID;
        cfg.dataset        = [files(nF).folder filesep files(nF).name];


        cfg.demean         = 'yes';
        cfg.lpfilter       = 'yes';        % enable low-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord      = 4;
        cfg.hpfreq         = 1; %!!! Filtering at 1Hz to improve the ICA decomposition
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50 100];     % set up the frequencies for notch filtering

        cfg.reref          = 'yes';
        cfg.refchannel     = 'all';

        data               = ft_preprocessing(cfg); % read raw data

        cfg.hpfreq         = 0.1; %
        data_eeg               = ft_preprocessing(cfg); % read raw data

        cfg=[];
        cfg.length  = 10; % Length of each epoch in seconds
        cfg.overlap = 0;  % No overlap between epochs
        data_epoched = ft_redefinetrial(cfg, data);
        data_eeg_epoched = ft_redefinetrial(cfg, data_eeg);

        %%% rename channels
        mylabels = data.label;
        for nCh = 1:length(mylabels)
            findspace = findstr(mylabels{nCh},' ');
            if isempty(findspace)
                newlabels{nCh} = mylabels{nCh};
            else
                if ismember(mylabels{nCh}(1),{'1','2','3','4','5','6','7','8','9'})
                    newlabels{nCh} = mylabels{nCh}(findspace+1:end);
                else
                    newlabels{nCh} = mylabels{nCh}(1:findspace-1);
                end
            end
        end
        cfg=[];
        cfg.channel          = data.label; %hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG'))) & (cellfun(@isempty,regexp(hdr.label,'EMG'))) & (cellfun(@isempty,regexp(hdr.label,'ECG'))) & (cellfun(@isempty,regexp(hdr.label,'Mic')))));
        cfg.montage.labelold = data.label;
        cfg.montage.labelnew = newlabels;
        cfg.montage.tra      = eye(length(data.label));
        data_epoched = ft_preprocessing(cfg, data_epoched); %causing an error because they're mutually exclusive
        data_eeg_epoched = ft_preprocessing(cfg, data_eeg_epoched); %causing an error because they're mutually exclusive

        % Automating removal of noisy channels
        std_vec = [];
        kurt_vec = [];
        for k = 1:length(data_eeg_epoched.trial)
            std_vec =[std_vec log(std(data_eeg_epoched.trial{k},[],2))];
            kurt_vec =[kurt_vec log(kurtosis(data_eeg_epoched.trial{k},[],2))];
        end
        all_std_vec = (reshape(std_vec,1,numel(std_vec)));
        all_kurt_vec = (reshape(kurt_vec,1,numel(kurt_vec)));
        all_std_vec_frontal = (reshape(std_vec(find_trials(newlabels,'F*'),:),1,numel(std_vec(find_trials(newlabels,'F*'),:))));
        all_kurt_vec_frontal = (reshape(kurt_vec(find_trials(newlabels,'F*'),:),1,numel(kurt_vec(find_trials(newlabels,'F*'),:))));
        all_std_vec_nonfront = (reshape(std_vec(setdiff(1:length(newlabels),find_trials(newlabels,'F*')),:),1,numel(std_vec(setdiff(1:length(newlabels),find_trials(newlabels,'F*')),:))));
        all_kurt_vec_nonfront = (reshape(kurt_vec(setdiff(1:length(newlabels),find_trials(newlabels,'F*')),:),1,numel(kurt_vec(setdiff(1:length(newlabels),find_trials(newlabels,'F*')),:))));

        thr_badElec_std=(mean(all_std_vec_nonfront)+3*std(all_std_vec_nonfront))*ones(length(newlabels),length(data_eeg_epoched.trial));
        thr_badElec_std(find_trials(newlabels,'F*'),:)=mean(all_std_vec_frontal)+3*std(all_std_vec_frontal);
        thr_badElec_kurt=(mean(all_kurt_vec_nonfront)+3*std(all_kurt_vec_nonfront))*ones(length(newlabels),length(data_eeg_epoched.trial));
        thr_badElec_kurt(find_trials(newlabels,'F*'),:)=mean(all_kurt_vec_frontal)+3*std(all_kurt_vec_frontal);
        badCh_std_mat = (std_vec>thr_badElec_std);
        badCh_kur_mat = (kurt_vec>thr_badElec_kurt);

        % we can define a bad channel as a channel with more than 50% of
        % probes above the threshold
        badCh_std=find(mean(badCh_std_mat')>0.5);
        badCh_kur=find(mean(badCh_kur_mat')>0.5);
        badCh=find(mean(((badCh_std_mat+badCh_kur_mat)~=0)')>0.5);
        % we can define a bad epoch as an epoch in which more than 20% of
        % channels are above threshold
        % probes above the threshold
        badTr_std=find(mean(badCh_std_mat)>0.2);
        badTr_kur=find(mean(badCh_kur_mat)>0.2);
        badTr=find(mean(((badCh_std_mat+badCh_kur_mat)~=0))>0.2);
        %
        %         if strcmp(SubID, 'C017')==1 %EP
        %             badCh = [badCh 17]; %This is forcing an interpolation of TP9 by adding it into the vector
        %         end
        badChannels=badCh;
        badTrials=badTr;

        nFc=nFc+1;
        badChannels_badTrials_info{nFc,1}=SubID;
        badChannels_badTrials_info{nFc,2}=badCh_std;
        badChannels_badTrials_info{nFc,3}=badCh_kur;
        badChannels_badTrials_info{nFc,4}=badChannels;
        badChannels_badTrials_info{nFc,5}=badTr_std;
        badChannels_badTrials_info{nFc,6}=badTr_kur;
        badChannels_badTrials_info{nFc,7}=badTrials;

        fprintf('Processing %s\n', SubID); % EP - debugging to check which participant is being processed
        %         disp(badChannels_badTrials_info(:, 1)); % EP - debugging to see if variable is being updated

        if ~isempty(badChannels)
            fprintf('... ... interpolating %g channels\n',length(badChannels))
            % find neighbours
            cfg=[];
            cfg.method        = 'triangulation';
            cfg.layout        = layout;
            cfg.feedback      = 'no';
            cfg.channel = layout.label;
            [neighbours] = ft_prepare_neighbours(cfg);

            % interpolate channels
            cfg=[];
            cfg.method         = 'weighted';
            cfg.badchannel     = layout.label(badChannels);
            cfg.missingchannel = [];
            cfg.neighbours     = neighbours;
            cfg.trials         = 'all';
            cfg.layout         = layout;
            cfg.channel = layout.label;
            [data_epoched] = ft_channelrepair(cfg, data_epoched); %crashing here
        end

        cfg=[];
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        data_epoched = ft_preprocessing(cfg,data_epoched);
        data_eeg_epoched = ft_preprocessing(cfg,data_eeg_epoched);

        if ~isempty(badTrials)
            cfg=[];
            cfg.trials    = setdiff(1:length(data_epoched.trial),badTrials);
            data_epoched = ft_redefinetrial(cfg,data_epoched);

            cfg=[];
            cfg.trials    = setdiff(1:length(data_eeg_epoched.trial),badTrials);
            data_eeg_epoched = ft_redefinetrial(cfg,data_eeg_epoched);
        end

        EEG = fieldtrip2eeglab(data_epoched);
        EEG = pop_chanedit(EEG, 'lookup', fullfile(path_eeglab, 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'));
        EEG = eeg_checkset(EEG); %This checks if current channel no. has the same amount as channel locations. If not, it deletes channel locations
        EEG_ica = pop_runica(EEG, 'icatype', 'runica'); %Runs ICA
        EEG_icalabels = pop_iclabel(EEG_ica,'default'); %Automates detectio of bad ICA components


        ICA_classification=EEG_icalabels.etc.ic_classification.ICLabel.classifications;
        ICA_classification=array2table(ICA_classification,'VariableNames',EEG_icalabels.etc.ic_classification.ICLabel.classes);

        ICA_classification.SubID=nan(size(ICA_classification,1),1);
        ICA_classification.SubID=categorical(ICA_classification.SubID);
        ICA_classification.SubID=repmat(SubID,size(ICA_classification,1),1);
        ICA_classification.Comp=(1:size(ICA_classification,1))';

        save([preproc_path filesep 'comp_i_RS_' SubID],'EEG_ica','EEG_icalabels','ICA_classification')
        all_ICA_classification=[all_ICA_classification ; ICA_classification]; %EP - moved this

        run('../MWMB_ADHD_elec_layout.m') % EP - moved this
%         figure('visible','off');
%         for nComp=1:16% size(EEG_ica.icawinv,2)
%             subplot(4,4,nComp)
%             simpleTopoPlot_ft(EEG_ica.icawinv(:,nComp), layout,'on',[],0,1); colorbar;
%             [maxValue, maxIndex] = max(table2array(ICA_classification(nComp, 1:7))); %EP
%             thisLabel = ICA_classification.Properties.VariableNames{maxIndex};
%             %thisLabel=ICA_classification.Properties.VariableNames(find(table2array(ICA_classification(nComp,1:7))==max(table2array(ICA_classification(nComp,1:7)))));
%             title(sprintf('%s: %1.3f', thisLabel, maxValue));
%             %title(sprintf('%s: %1.3f',thisLabel{1},max(table2array(ICA_classification(nComp,:)))));
%         end
%         savefig(gcf,[preproc_path filesep 'comp_i_RS_' SubID '.fig'])
%         close(gcf)

        data_eeg_epoched.hdr.Fs=data_eeg_epoched.fsample;
        data_eeg_epoched.hdr.nChans=length(data_eeg_epoched.label);
        data_eeg_epoched.hdr.label=data_eeg_epoched.label;
        data_eeg_epoched.hdr.nSamplesPre=data_eeg_epoched.time{1}(1)*data_eeg_epoched.fsample;
        data_eeg_epoched.hdr.nSamples=length(data_eeg_epoched.time{1})*length(data_eeg_epoched.time);
        data_eeg_epoched.hdr.nTrials=length(data_eeg_epoched.time);
        eeglab_eeg_data = fieldtrip2eeglab(data_eeg_epoched);

        rejected_comps = ICA_classification.Comp(ICA_classification.Eye>0.9 | ICA_classification.Heart>0.8); % reject eye component with proba over 0.95 and heart over 0.8
        fprintf('... ... %g bad components rejected\n',length(rejected_comps))
        EEG_ica.data = eeglab_eeg_data.data;
        EEG_clean = pop_subcomp(EEG_ica, rejected_comps);
        EEG_clean = eeg_checkset(EEG_clean);
        badCompo=ICA_classification(rejected_comps,:);
        %                     figure;
        %                     plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'Fp1'),:,1))','r');
        %                     hold on
        %                     plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'Fp2'),:,1))','b');
        %                     plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Fp1'),:,1))','r--');
        %                     plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'Fp2'),:,1))','b--');
        %
        %                         figure;
        %                     plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'O1'),:,1))','r');
        %                     hold on
        %                     plot(squeeze(EEG_ica.data(match_str({EEG_ica.chanlocs.labels},'O2'),:,1))','b');
        %                     plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'O1'),:,1))','r--');
        %                     plot(squeeze(EEG_clean.data(match_str({EEG_ica.chanlocs.labels},'O2'),:,1))','b--');
        ori_data=data_eeg_epoched;
        data = eeglab2fieldtrip(EEG_clean, 'raw');
        %         figure;
        %         plot(squeeze(data.trial{2}(match_str(data.label,'Fz'),:)),'r');
        %         hold on
        %         plot(squeeze(data.trial{2}(match_str(data.label,'Cz'),:)),'b');
        %         plot(squeeze(ori_data.trial{2}(match_str(data.label,'Fz'),:)),'r--');
        %         hold on
        %         plot(squeeze(ori_data.trial{2}(match_str(data.label,'Cz'),:)),'b--');
        %
        save([preproc_path filesep 'clean_i_RS_' SubID '.mat'],'data','badCompo');
        all_badCompo=[all_badCompo ; badCompo];


        %%% Get power
        cfg                             = [];
        cfg.foilim                      = [2 40]; %[1 30];
        cfg.pad                         = 'nextpow2';
        cfg.tapsmofrq                   = 0.5;
        cfg.method                      = 'mtmfft';
        cfg.output                      = 'fooof_aperiodic';
        cfg.fooof.max_peaks             = 4;
        cfg.fooof.proximity_threshold   = 1;

        cfg.keeptrials                  = 'no';
        cfg.output                      = 'pow';
        pow = ft_freqanalysis(cfg, data);
        pow_noICA = ft_freqanalysis(cfg, data_eeg_epoched);
        save([preproc_path filesep 'clean_pow_RS_' SubID '.mat'],'pow','pow_noICA');
    end
    toc;
end



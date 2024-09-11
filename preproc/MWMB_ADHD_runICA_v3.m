%% Changes from v2 - flagging of badChannels (changed from 4std to 3std) and therefore saved files all have ver2 in them
%% 04/09/24 - EP: Using this script to debug/investigate why where probes are going missing but NOTE: v2 Thomas fixed so that there's a different threshold for frontal electrodes 
%% NOT USING THIS SCRIPT FOR RUNNING ICAs!!
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

if exist([preproc_path filesep 'all_badChannels_badProbes_v2.mat'])%EP - load bad channels table 
    load([preproc_path filesep 'all_badChannels_badProbes_v2.mat']);
    startIndex = size(badChannels_badTrials_info_v2, 1) + 1;
else
    badChannels_badTrials_info_v2 = {};
    startIndex = 1;
end 

%% Loop across files
RS = ["R1", "R2"];

redo=1;
nFc = startIndex-1;
all_ICA_classification=[];
for nF=1:length(eeg_files)
    if startsWith(eeg_files(nF).name, '._') % EP - Skip this file if it starts with dot underline.
        continue; %  EP - Jump to the bottom of the loop.
    end

    if contains(eeg_files(nF).name,RS) %To skip resting state files
        continue
    end
%         if ~contains(eeg_files(nF).name,'ID-C038') %To redo specific a subject
%         continue
%     end

    %%% load the data
    SubInfo=split(eeg_files(nF).name,'-');
    SubID=SubInfo{2}(1:end-4);

        if ~contains(SubID,'A039') %To debug - looking for participants who had probes where no SW were detected
        continue
    end

    if redo==1 || exist([preproc_path filesep 'comp2_i_probe_' SubID '.mat'])==0 || ~any(strcmp(badChannels_badTrials_info_v2(:,1), SubID))% To skip already preprocessed files
        fprintf('... working on %s\n',[eeg_files(nF).name])

        %%% minimal preprocessing
        cfg=[];
        cfg.SubID          = SubID;
        cfg.dataset        = [eeg_files(nF).folder filesep eeg_files(nF).name];


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

        cfg.trialfun            = 'MWMB_ADHD_probefun';
        %         cfg.table               = table;
        cfg.SubID               = SubID;
        cfg.dataset             = [eeg_files(nF).folder filesep eeg_files(nF).name];
        cfg.trialdef.prestim    = 25;
        cfg.trialdef.poststim   = 2;
        cfg = ft_definetrial(cfg);

        data               = ft_preprocessing(cfg); % read raw data

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
        data = ft_preprocessing(cfg, data); %causing an error because they're mutually exclusive

        % Automating removal of noisy channels
        std_vec = [];
        kurt_vec = [];
        for k = 1:length(data.trial)
            std_vec =[std_vec log(std(data.trial{k},[],2))];
            kurt_vec =[kurt_vec log(kurtosis(data.trial{k},[],2))];
        end
        all_std_vec = (reshape(std_vec,1,numel(std_vec)));
        all_kurt_vec = (reshape(kurt_vec,1,numel(kurt_vec)));
        badCh_std_mat = (std_vec>(mean(all_std_vec)+3*std(all_std_vec))); % in version 2, ths was 4* std
        badCh_kur_mat = (kurt_vec>(mean(all_kurt_vec)+3*std(all_kurt_vec))); % in version 2, ths was 4* std

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
        if strcmp(SubID, 'C017')==1 %EP
            badCh = [badCh 17]; %This is forcing an interpolation of TP9 by adding it into the vector
        end
        badChannels=badCh;
        badTrials=badTr;

        nFc=nFc+1;
        badChannels_badTrials_info_v2{nFc,1}=SubID; % need to edit to load this info up if it has been processed already
        badChannels_badTrials_info_v2{nFc,2}=badCh_std;
        badChannels_badTrials_info_v2{nFc,3}=badCh_kur;
        badChannels_badTrials_info_v2{nFc,4}=badChannels;
        badChannels_badTrials_info_v2{nFc,5}=badTr_std;
        badChannels_badTrials_info_v2{nFc,6}=badTr_kur;
        badChannels_badTrials_info_v2{nFc,7}=badTrials;

        fprintf('Processing %s\n', SubID); % EP - debugging to check which participant is being processed
        disp(badChannels_badTrials_info_v2(:, 1));; % EP - debugging to see if variable is being updated
        
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
            [data] = ft_channelrepair(cfg, data); 
        end

        cfg=[];
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        data = ft_preprocessing(cfg,data);

        if ~isempty(badTrials)
            cfg=[];
            cfg.trials    = setdiff(1:length(data.trial),badTrials);
            data = ft_redefinetrial(cfg,data);
        end

        EEG = fieldtrip2eeglab(data);
        EEG = pop_chanedit(EEG, 'lookup', fullfile(path_eeglab, 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'));
        EEG = eeg_checkset(EEG); %This checks if current channel no. has the same amount as channel locations. If not, it deletes channel locations
        EEG_ica = pop_runica(EEG, 'icatype', 'runica'); %Runs ICA
        EEG_icalabels = pop_iclabel(EEG_ica,'default'); %Automates detection of bad ICA components


        ICA_classification=EEG_icalabels.etc.ic_classification.ICLabel.classifications;
        ICA_classification=array2table(ICA_classification,'VariableNames',EEG_icalabels.etc.ic_classification.ICLabel.classes);

        ICA_classification.SubID=nan(size(ICA_classification,1),1);
        ICA_classification.SubID=categorical(ICA_classification.SubID);
        ICA_classification.SubID=repmat(SubID,size(ICA_classification,1),1);
        ICA_classification.Comp=(1:size(ICA_classification,1))';

        save([preproc_path filesep 'comp2_i_probe_' SubID],'EEG_ica','EEG_icalabels','ICA_classification')
        all_ICA_classification=[all_ICA_classification ; ICA_classification]; %EP - moved this

        run('../MWMB_ADHD_elec_layout.m') % EP - moved this
        figure('visible','off');
        for nComp=1:16% size(EEG_ica.icawinv,2)
            subplot(4,4,nComp)
            simpleTopoPlot_ft(EEG_ica.icawinv(:,nComp), layout,'on',[],0,1); colorbar;
            [maxValue, maxIndex] = max(table2array(ICA_classification(nComp, 1:7))); %EP
            thisLabel = ICA_classification.Properties.VariableNames{maxIndex};
            %thisLabel=ICA_classification.Properties.VariableNames(find(table2array(ICA_classification(nComp,1:7))==max(table2array(ICA_classification(nComp,1:7)))));
            title(sprintf('%s: %1.3f', thisLabel, maxValue));
            %title(sprintf('%s: %1.3f',thisLabel{1},max(table2array(ICA_classification(nComp,:)))));
        end
        savefig(gcf,[preproc_path filesep 'comp2_i_probe_' SubID '.fig'])
        close(gcf)

    else
        load([preproc_path filesep 'comp2_i_probe_' SubID]) %EP
        fprintf('... working on %s\n',[eeg_files(nF).name]) %EP
        all_ICA_classification=[all_ICA_classification ; ICA_classification]; %EP

        run('../MWMB_ADHD_elec_layout.m') % EP - moved this
        figure('visible','off');
        for nComp=1:16% size(EEG_ica.icawinv,2)
            subplot(4,4,nComp)
            simpleTopoPlot_ft(EEG_ica.icawinv(:,nComp), layout,'on',[],0,1); colorbar;
            [maxValue, maxIndex] = max(table2array(ICA_classification(nComp, 1:7))); %EP
            thisLabel = ICA_classification.Properties.VariableNames{maxIndex};
            %thisLabel=ICA_classification.Properties.VariableNames(find(table2array(ICA_classification(nComp,1:7))==max(table2array(ICA_classification(nComp,1:7)))));
            title(sprintf('%s: %1.3f', thisLabel, maxValue));
            %title(sprintf('%s: %1.3f',thisLabel{1},max(table2array(ICA_classification(nComp,:)))));
        end
        savefig(gcf,[preproc_path filesep 'comp2_i_probe_' SubID '.fig'])
        close(gcf)

    end
end

save([preproc_path filesep 'all_badChannels_badProbes_v2'],'badChannels_badTrials_info_v2')

% % EP - Commented out 238-250 to move them into the loop
% run('../MWMB_ADHD_elec_layout.m')
% figure('visible','off');
% for nComp=1:16% size(EEG_ica.icawinv,2)
%     subplot(4,4,nComp)
%     simpleTopoPlot_ft(EEG_ica.icawinv(:,nComp), layout,'on',[],0,1); colorbar;
%     [maxValue, maxIndex] = max(table2array(ICA_classification(nComp, 1:7))); %EP
%     thisLabel = ICA_classification.Properties.VariableNames{maxIndex};
%     %thisLabel=ICA_classification.Properties.VariableNames(find(table2array(ICA_classification(nComp,1:7))==max(table2array(ICA_classification(nComp,1:7)))));
%     title(sprintf('%s: %1.3f', thisLabel, maxValue));
%     %title(sprintf('%s: %1.3f',thisLabel{1},max(table2array(ICA_classification(nComp,:)))));
% end
% savefig(gcf,[preproc_path filesep 'comp_i_probe_' SubID '.fig'])
% close(gcf)


%         rejected_comps = find(EEG.reject.gcompreject > 0);
%         EEG = pop_subcomp(EEG, rejected_comps);
%         EEG = eeg_checkset(EEG);

%         % convert back to Fieldtrip
%         curPath = pwd;
%         p = fileparts(which('ft_read_header'));
%         cd(fullfile(p, 'private'));
%         hdr = read_eeglabheader( EEG );
%         data = read_eeglabdata( EEG, 'header', hdr );
%         event = read_eeglabevent( EEG, 'header', hdr );
%
%         OUTEEG = pop_subcomp( EEG, components, plotag);
%
%         EEG_icalabels = pop_iclabel(EEG_ica,'default');
%
%         %defining epochs - moved this after ICA because the script wasn't
%         %running due to diff trial numbers per block
%         cfg                     = [];
%         cfg.trialfun            = 'MWMB_ADHD_blockfun';
%         cfg.trialdef.prestim    = 1;
%         cfg.trialdef.poststim   = 1;
%         cfg                     = ft_definetrial(cfg);
%
%         save([preproc_path filesep 'Icfe_MWADHD_' SubID(1:end-4) '.mat'],'data','comp','rankICA','badChannels');

set(gcf,'Visible','on')
writetable(all_ICA_classification,[preproc_path filesep 'ICA_classification_allSubs_v2.csv'])
figure;
histogram(all_ICA_classification.Eye,0:0.05:1)
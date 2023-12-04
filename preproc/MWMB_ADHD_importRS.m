%%
clear all;
close all;

%%
if isempty(findstr(pwd,'thandrillon'))==0
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_data='/Users/thandrillon/Data/ADHD_MW/';
else
    
end

run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools))

files=dir([path_data filesep 'EEG' filesep '*_R*.eeg']);

%% loop on subjects
redo=0; % 1 to force re-import, 0 otherwise
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    if isempty(findstr(SubID,'ID-'))==0
        SubID(1:findstr(SubID,'ID-')+3)=[];
    end
    if isempty(findstr(SubID,'_R'))==0
        SubID(findstr(SubID,'_R'):end)=[];
    end
    tic;
    %     if ~strcmp(SubID,'A084') %%re-add this if you need to clip the data
    %         continue;
    %     end
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    if redo==1 || exist([path_data filesep 'Preproc' filesep 'fe_ft_RS1_' SubID '.mat'])==0
        
        hdr=ft_read_header([folder_name filesep file_name]);
        evt=ft_read_event([folder_name filesep file_name]);
        
        %%% Define RS1
        evt(find(cellfun(@isempty,{evt.value})==1)) = []; %% Elaine - added to remove empty cells in evt.value
        
        cfg=[];
        cfg.SubID               = SubID;
        cfg.dataset             = [folder_name filesep file_name];
        
        cfg.channel        = hdr.label(match_str(hdr.chantype,'eeg'));
        cfg.demean         = 'yes';
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord      = 4;
        cfg.hpfreq         = 0.1;
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering
        
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        
        dat                = ft_preprocessing(cfg); % read raw data
        
        cfgbs=[];
        cfgbs.resamplefs      = 250;
        cfgbs.detrend         = 'no';
        cfgbs.demean          = 'yes';
        cfgbs.baselinewindow  = [-0.2 0];
        data                  = ft_resampledata(cfgbs,dat); % read raw data
        save([path_data filesep 'Preproc' filesep 'fe_ft_' SubID],'data');
        
    end
    toc;
end



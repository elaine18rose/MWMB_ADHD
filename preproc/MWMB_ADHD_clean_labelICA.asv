clear all;
close all;

%% Paths
if isempty(findstr(pwd,'thandrillon'))==0
    path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    data_path='/Users/thandrillon/Data/Nir_NeuroMod/raw_data/';
    preproc_path='/Users/thandrillon/Data/Nir_NeuroMod/preproc/';
    path_detectSW='/Users/thandrillon/Data/Nir_NeuroMod/SW_detection/';
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

addpath((path_fieldtrip));
addpath(genpath(path_LSCPtools));
ft_defaults;
addpath(path_eeglab);
eeglab;

%% List all the available files
eeg_files=dir([data_path filesep '*.eeg']);

% Marie's: dir_probe=dir([path_preproc filesep 'probe_*.mat']);

%%
for nF=1:length(eeg_files)
    disp([eeg_files(nF).name]);
    % load data
    try
        load([preproc_path filesep eeg_files(nF).name]);
    catch
        Problem_Files{end+1}=eeg_files(nF).name;
        warning(sprintf('This file cannot be loaded',eeg_files(nF).name))
        continue;
    end
    
    % load([path_preproc filesep 'i_' dir_probe(nF).name])
    load([preproc_path filesep 'Icfe_MWADHD_' eeg_files(nF).name])
    
    cfg=[];
    cfg.reref      = 'yes';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg,data);
    
    cfg=[];
    cfg.toilim    = [-30 0];
    data = ft_redefinetrial(cfg,data);
    
    EEG = fieldtrip2eeglab(data);
    eloc = readlocs('/Users/elaine/Desktop/MATLAB_Functions/fieldtrip/template/layout/acticap-64ch-standard2.mat'); %%eloc = readlocs('chanlocs.ced'); % Channel location - right now it's 28 channels when we need 64
    toRemove = false(1, length(eloc));
    for nCh = 1:length(eloc)
        if contains(eloc(nCh).labels,'Gnd')|| contains(eloc(nCh).labels,'Ref')
            toRemove(nCh) = true;
        end
    end
    eloc(toRemove) = []; %Remove Grnd and Ref electrodes

    EEG.chanlocs=eloc;
    EEG = eeg_checkset(EEG);

    %EEG = pop_runica(EEG, 'icatype', 'runica'); % this was in runICA script
    EEG.icawinv=comp.topo;
    EEG.icaweights=comp.unmixing;
    
%EEG = pop_iclabel(EEG, 'default'); % see function help message this was
%also in run ICA script 
rejected_comps = find(EEG.reject.gcompreject > 0);
EEG = pop_subcomp(EEG, rejected_comps);
EEG = eeg_checkset(EEG);

% convert back to Fieldtrip
curPath = pwd;
p = fileparts(which('ft_read_header'));
cd(fullfile(p, 'private'));
hdr = read_eeglabheader( EEG );
data = read_eeglabdata( EEG, 'header', hdr );
event = read_eeglabevent( EEG, 'header', hdr );

OUTEEG = pop_subcomp( EEG, components, plotag);


EEG_icalabels = pop_iclabel(EEG_ica,'default');

end


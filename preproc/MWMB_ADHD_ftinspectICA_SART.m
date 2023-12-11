%%
if isempty(findstr(pwd,'thandrillon'))==0
    path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    data_path='/Users/thandrillon/Data/Nir_NeuroMod/raw_data/';
    preproc_path='/Users/thandrillon/Data/Nir_NeuroMod/preproc/';
    path_detectSW='/Users/thandrillon/Data/Nir_NeuroMod/SW_detection/';
else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/desktop/MATLAB_Functions/fieldtrip/';
    data_path = '/Volumes/Seagate/MWMB_ADHD_SART/EEG/';
    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    path_detectSW = '/Volumes/Seagate/MWMB_ADHD_SART/SW_detection/';

    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(path_fieldtrip);
ft_defaults;

%% choose and load subject
List_Subj=dir([preproc_path filesep 'Icfe_ft_MWADHD_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
load([preproc_path filesep ListNames{pick}])
oridata=data;

%% Layout
cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.channel=data.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

%% Plot Components
figure;
cfg = [];
cfg.component = 1:32;       % specify the component(s) that should be plotted
cfg.layout    = layout; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp)
set(gcf,'Position',[1           1        1871         984]);

figure;
cfg = [];
cfg.component = 33:length(comp.label);       % specify the component(s) that should be plotted
cfg.layout    = layout;                      % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp)
set(gcf,'Position',[1           1        1871         984]);


%% reject trials
pickComponents=input('Select component you want to plot: ');

%%
cfg          = [];
cfg.channel  = pickComponents; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = layout; % specify the layout file that should be used for plotting
cfg.allowoverlap='true';
cfg.continuous='no';
ft_databrowser(cfg, comp);
fprintf('... working on %s\n',ListNames{pick})


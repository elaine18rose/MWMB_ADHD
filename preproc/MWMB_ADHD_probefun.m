function trl = MWMB_ADHD_probefun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
trl = [];

evt=ft_read_event(cfg.dataset);
%               trig_start          =1;  %S
%             trig_end            =4;  %E
%             trig_startBlock     =5;  %B
%             trig_endBlock       =8;  %K
%             trig_startGO        =9;  %T
%             trig_startNOGO      =10; %T
%             trig_startQuestion  =12; %Q
%             trig_probestart     =13; %P
%             trig_probeend       =16; %C
evt(find(cellfun(@isempty,{evt.value})==1)) = []; %% Elaine - added to remove empty cells in evt.value
evt(find_trials({evt.type},'New Segment')) = []; %% Elaine - added to remove empty cells in evt.value
startProbeIdx=find_trials({evt.value},'S 13');
if strcmp(cfg.SubID,'A015')
    startProbeIdx=[startProbeIdx 1012]; % wrong triggers
    startProbeIdx=sort(startProbeIdx);
elseif strcmp(cfg.SubID,'A050')
    startProbeIdx=[startProbeIdx 1207]; % probe occurred during a buffer overload
    startProbeIdx=sort(startProbeIdx);
elseif strcmp(cfg.SubID,'C042')
    startProbeIdx=[startProbeIdx 2009]; % probe occurred during a buffer overload
    startProbeIdx=sort(startProbeIdx);
elseif strcmp(cfg.SubID,'A038')
    startProbeIdx=[startProbeIdx 377 443 486]; % probe occurred during a buffer overload
    startProbeIdx=sort(startProbeIdx);
elseif strcmp(cfg.SubID,'A039')
    startProbeIdx=[startProbeIdx 428 481 529 587 647 ...
        693 740 799 854 916 976 1083 1138 1184 1236 ...
        1301 1366 1411 1478 1534 1585 1641 1689 1794 1860 ...
        1921 1975 2036 2099 2167 2235 2343 2391 2436 2490];
    startProbeIdx=sort(startProbeIdx); % EP - wrong triggers; start triggers became S 8s
    %%%Code below is to help with troubelshooting ppts with diff number of probes and triggers in EEG and behav file
    % load('/Users/thandrillon/Data/ADHD_MW/Behaviour/wanderIM_behavres_A039_29May2024-1210.mat')
%     diff_probeend=diff(probe_res(:,3));
%     evt_samples=[evt.sample];
%     [diff(evt_samples(startProbeIdx))/hdr.Fs ; diff_probeend']
end
evt_samples=[evt.sample];
for i=1:length(startProbeIdx)
    % add this to the trl definition
    begsample     = evt_samples(startProbeIdx(i))-cfg.trialdef.prestim*hdr.Fs;
    endsample     = evt_samples(startProbeIdx(i))+cfg.trialdef.poststim*hdr.Fs;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end

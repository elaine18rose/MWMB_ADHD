function [trl,myevt] = MWMB_ADHD_trialfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
trl = [];
myevt=[];
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
startTrialIdx=unique([find_trials({evt.value},{'S  9'}) find_trials({evt.value},{'S 10'})]);
myevt=evt(startTrialIdx);
evt_samples=[evt.sample];
for i=1:length(startTrialIdx)
    % add this to the trl definition
    begsample     = evt_samples(startTrialIdx(i))-cfg.trialdef.prestim*hdr.Fs;
    endsample     = evt_samples(startTrialIdx(i))+cfg.trialdef.poststim*hdr.Fs;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end

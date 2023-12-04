function trl = MWMB_ADHD_blockfun(cfg)

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
startBlockIdx=find_trials({evt.value},'S  5');
endBlockIdx=find_trials({evt.value},'S  8');

if length(startBlockIdx)~=length(endBlockIdx)
    return;
end
evt_samples=[evt.sample];
for i=1:length(startBlockIdx)
    % add this to the trl definition
    begsample     = evt_samples(startBlockIdx(i));
    endsample     = evt_samples(endBlockIdx(i));
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end

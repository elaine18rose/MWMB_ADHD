% FileName	RS	Bad_Channels
Labels={'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10','P7','P3','Pz','P4','P8','PO9','O1','Oz','O2','PO10','AF7','AF3','AF4','AF8','F5','F1','F2','F6','FT9','FT7','FC3','FC4','FT8','FT10','C5','C1','C2','C6','TP7','CP3','CPz','CP4','TP8','P5','P1','P2','P6','PO7','PO3','POz','PO4','PO8'};

cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.channel = Labels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
orderchan=[];
for nCh=1:length(Labels)
    orderchan(nCh)=match_str(layout.label,Labels{nCh});
end
orderchan(length(Labels)+1)=length(Labels)+1;
orderchan(length(Labels)+2)=length(Labels)+2;
layout.label=layout.label(orderchan);
layout.pos=layout.pos(orderchan,:);
layout.width=layout.width(orderchan);
layout.height=layout.height(orderchan);

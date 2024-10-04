function [all_clus]=get_clusterperm_lme_bygroup(model_est,clus_alpha,montecarlo_alpha,totperm,neighbours,limSize)

if nargin<6 
    limSize=[]; % ???
end
if size(model_est{2},2)==6 % ???
    for nCond=1:3 % For each drug condition (being compared to placebo)
        real_clus=cell(1,2);
        for nsign=1:2 % For positive and negative relationships
            if nsign==1
                sig_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)>0 & model_est{1}(:,5)==nCond,1); % Elec names for significant elecs with positive effects and group of interest
                tstat_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)>0 & model_est{1}(:,5)==nCond,3); % t values as above
            else
                sig_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)<0 & model_est{1}(:,5)==nCond,1); % As above but negative
                tstat_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)<0 & model_est{1}(:,5)==nCond,3); % As above but negative
            end
            if ~isempty(sig_elec) % If there are significant electrodes
                real_clus{nsign}{1,1}={neighbours(sig_elec(1)).label}; % Get the names of its neighbours
                real_clus{nsign}{1,3}=tstat_elec(1); % And its t value
                real_clus{nsign}{1,2}=neighbours(sig_elec(1)).neighblabel'; % And its neighbours? (how is this different?)
                for nEl=2:length(sig_elec) % For all significant electrodes
                    clus_found=0;
                    for nclus=1:size(real_clus{nsign},1) % For all clusters
                        if sum(ismember(real_clus{nsign}{nclus,1},neighbours(sig_elec(nEl)).neighblabel'))~=0 % If 
                            real_clus{nsign}{nclus,1}=unique([real_clus{nsign}{nclus,1} ; {neighbours(sig_elec(nEl)).label}]);
                            real_clus{nsign}{nclus,3}=[real_clus{nsign}{nclus,3} ; tstat_elec(nEl)];
                            real_clus{nsign}{nclus,2}=unique([real_clus{nsign}{nclus,2} ; neighbours(sig_elec(nEl)).neighblabel']);
                            clus_found=1;
                            break;
                        end
                    end
                    if clus_found==0
                        num_clus=size(real_clus{nsign},1);
                        real_clus{nsign}{num_clus+1,1}= {neighbours(sig_elec(nEl)).label};
                        real_clus{nsign}{num_clus+1,3}=tstat_elec(nEl);
                        real_clus{nsign}{num_clus+1,2}=neighbours(sig_elec(nEl)).neighblabel';
                    end
                end
                if ~isempty(limSize)
                    clus_discard=[];
                    for m=1:size(real_clus{nsign},1)
                        if size(real_clus{nsign}{m,1},1)<=limSize
                            clus_discard=[clus_discard m];
                        end
                    end
                    real_clus{nsign}(clus_discard,:)=[];
                end
                clus_merged=[];
                for nclus=1:size(real_clus{nsign},1)
                    if ismember(nclus,clus_merged)
                        continue;
                    end
                    for nclus2=nclus+1:size(real_clus{nsign},1)
                        if sum(ismember(real_clus{nsign}{nclus2,1},real_clus{nsign}{nclus,2}))~=0
                            real_clus{nsign}{nclus,3}=[real_clus{nsign}{nclus,3} ; real_clus{nsign}{nclus2,3}(~ismember(real_clus{nsign}{nclus2,1},real_clus{nsign}{nclus,1}))];
                            real_clus{nsign}{nclus,1}=unique([real_clus{nsign}{nclus,1} ; real_clus{nsign}{nclus2,1}]);
                            real_clus{nsign}{nclus,2}=unique([real_clus{nsign}{nclus,2} ; real_clus{nsign}{nclus2,2}]);
                            clus_merged=[clus_merged nclus2];
                        end
                    end
                end
                real_clus{nsign}(clus_merged,:)=[];
            end
        end
        
        perm_clus=cell(2,totperm);
        for nperm=1:totperm
            for nsign=1:2
                if nsign==1
                    sig_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)>0 & model_est{2}(:,5)==nperm & model_est{2}(:,6)==nCond,1);
                    tstat_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)>0 & model_est{2}(:,5)==nperm & model_est{2}(:,6)==nCond,3);
                else
                    sig_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)<0 & model_est{2}(:,5)==nperm & model_est{2}(:,6)==nCond,1);
                    tstat_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)<0 & model_est{2}(:,5)==nperm & model_est{2}(:,6)==nCond,3);
                end
                if ~isempty(sig_elec)
                    perm_clus{nsign,nperm}{1,1}={neighbours(sig_elec(1)).label};
                    perm_clus{nsign,nperm}{1,3}=tstat_elec(1);
                    perm_clus{nsign,nperm}{1,2}=neighbours(sig_elec(1)).neighblabel';
                    for nEl=2:length(sig_elec)
                        clus_found=0;
                        for nclus=1:size(perm_clus{nsign,nperm},1)
                            if sum(ismember(perm_clus{nsign,nperm}{nclus,1},neighbours(sig_elec(nEl)).neighblabel'))~=0
                                perm_clus{nsign,nperm}{nclus,1}=unique([perm_clus{nsign,nperm}{nclus,1} ; {neighbours(sig_elec(nEl)).label}]);
                                perm_clus{nsign,nperm}{nclus,3}=[perm_clus{nsign,nperm}{nclus,3} ; tstat_elec(nEl)];
                                perm_clus{nsign,nperm}{nclus,2}=unique([perm_clus{nsign,nperm}{nclus,2} ; neighbours(sig_elec(nEl)).neighblabel']);
                                clus_found=1;
                                break;
                            end
                        end
                        if clus_found==0
                            num_clus=size(perm_clus{nsign,nperm},1);
                            perm_clus{nsign,nperm}{num_clus+1,1}= {neighbours(sig_elec(nEl)).label};
                            perm_clus{nsign,nperm}{num_clus+1,3}=tstat_elec(nEl);
                            perm_clus{nsign,nperm}{num_clus+1,2}=neighbours(sig_elec(nEl)).neighblabel';
                        end
                    end
                    if ~isempty(limSize)
                        clus_discard=[];
                        for m=1:size(perm_clus{nsign,nperm},1)
                            if size(perm_clus{nsign,nperm}{m,1},1)<=limSize
                                clus_discard=[clus_discard m];
                            end
                        end
                        perm_clus{nsign,nperm}(clus_discard,:)=[];
                    end
                    clus_merged=[];
                    for nclus=1:size(perm_clus{nsign,nperm},1)
                        if ismember(nclus,clus_merged)
                            continue;
                        end
                        for nclus2=nclus+1:size(perm_clus{nsign,nperm},1)
                            if sum(ismember(perm_clus{nsign,nperm}{nclus2,1},perm_clus{nsign,nperm}{nclus,2}))~=0
                                perm_clus{nsign,nperm}{nclus,3}=[perm_clus{nsign,nperm}{nclus,3} ; perm_clus{nsign,nperm}{nclus2,3}(~ismember(perm_clus{nsign,nperm}{nclus2,1},perm_clus{nsign,nperm}{nclus,1}))];
                                perm_clus{nsign,nperm}{nclus,1}=unique([perm_clus{nsign,nperm}{nclus,1} ; perm_clus{nsign,nperm}{nclus2,1}]);
                                perm_clus{nsign,nperm}{nclus,2}=unique([perm_clus{nsign,nperm}{nclus,2} ; perm_clus{nsign,nperm}{nclus2,2}]);
                                clus_merged=[clus_merged nclus2];
                            end
                        end
                    end
                    perm_clus{nsign,nperm}(clus_merged,:)=[];
                end
            end
        end
        
        perm_tval=cell(1,2);
        for nsign=1:2
            for nperm=1:totperm
                temp_tval=[];
                for nclus=1:size(perm_clus{nsign,nperm},1)
                    temp_tval=[temp_tval sum(perm_clus{nsign,nperm}{nclus,3})];
                end
                if nsign==1
                    if isempty(temp_tval)
                        perm_tval{nsign}=[perm_tval{nsign} 0];
                    else
                        perm_tval{nsign}=[perm_tval{nsign} max(temp_tval)];
                    end
                else
                    if isempty(temp_tval)
                        perm_tval{nsign}=[perm_tval{nsign} 0];
                    else
                        perm_tval{nsign}=[perm_tval{nsign} min(temp_tval)];
                    end
                end
            end
        end
        
        all_clus{nCond}=[]; nc=0;
        for nsign=1:2
            for nclus=1:size(real_clus{nsign},1)
                this_tval=sum(real_clus{nsign}{nclus,3});
                if nsign==1
                    if mean(this_tval<perm_tval{nsign})<montecarlo_alpha
                        nc=nc+1;
                        all_clus{nCond}{nc}={'pos' real_clus{nsign}{nclus,1} this_tval sum(this_tval<perm_tval{nsign})/totperm};
                    end
                else
                    if mean(this_tval>perm_tval{nsign})<montecarlo_alpha
                        nc=nc+1;
                        all_clus{nCond}{nc}={'neg' real_clus{nsign}{nclus,1} this_tval sum(this_tval>perm_tval{nsign})/totperm};
                    end
                end
            end
        end
    end
else
    real_clus=cell(1,2);
    for nsign=1:2
        if nsign==1
            sig_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)>0 ,1);
            tstat_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)>0 ,3);
        else
            sig_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)<0 ,1);
            tstat_elec=model_est{1}(model_est{1}(:,4)<clus_alpha & model_est{1}(:,3)<0 ,3);
        end
        if ~isempty(sig_elec)
            real_clus{nsign}{1,1}={neighbours(sig_elec(1)).label};
            real_clus{nsign}{1,3}=tstat_elec(1);
            real_clus{nsign}{1,2}=neighbours(sig_elec(1)).neighblabel';
            for nEl=2:length(sig_elec)
                clus_found=0;
                for nclus=1:size(real_clus{nsign},1)
                    if sum(ismember(real_clus{nsign}{nclus,1},neighbours(sig_elec(nEl)).neighblabel'))~=0
                        real_clus{nsign}{nclus,1}=unique([real_clus{nsign}{nclus,1} ; {neighbours(sig_elec(nEl)).label}]);
                        real_clus{nsign}{nclus,3}=[real_clus{nsign}{nclus,3} ; tstat_elec(nEl)];
                        real_clus{nsign}{nclus,2}=unique([real_clus{nsign}{nclus,2} ; neighbours(sig_elec(nEl)).neighblabel']);
                        clus_found=1;
                        break;
                    end
                end
                if clus_found==0
                    num_clus=size(real_clus{nsign},1);
                    real_clus{nsign}{num_clus+1,1}= {neighbours(sig_elec(nEl)).label};
                    real_clus{nsign}{num_clus+1,3}=tstat_elec(nEl);
                    real_clus{nsign}{num_clus+1,2}=neighbours(sig_elec(nEl)).neighblabel';
                end
            end
            if ~isempty(limSize)
                clus_discard=[];
                for m=1:size(real_clus{nsign},1)
                    if size(real_clus{nsign}{m,1},1)<=limSize
                        clus_discard=[clus_discard m];
                    end
                end
                real_clus{nsign}(clus_discard,:)=[];
            end
            clus_merged=[];
            for nclus=1:size(real_clus{nsign},1)
                if ismember(nclus,clus_merged)
                    continue;
                end
                for nclus2=nclus+1:size(real_clus{nsign},1)
                    if sum(ismember(real_clus{nsign}{nclus2,1},real_clus{nsign}{nclus,2}))~=0
                        real_clus{nsign}{nclus,3}=[real_clus{nsign}{nclus,3} ; real_clus{nsign}{nclus2,3}(~ismember(real_clus{nsign}{nclus2,1},real_clus{nsign}{nclus,1}))];
                        real_clus{nsign}{nclus,1}=unique([real_clus{nsign}{nclus,1} ; real_clus{nsign}{nclus2,1}]);
                        real_clus{nsign}{nclus,2}=unique([real_clus{nsign}{nclus,2} ; real_clus{nsign}{nclus2,2}]);
                        clus_merged=[clus_merged nclus2];
                    end
                end
            end
            real_clus{nsign}(clus_merged,:)=[];
        end
    end
    
    perm_clus=cell(2,totperm);
    for nperm=1:totperm
        for nsign=1:2
            if nsign==1
                sig_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)>0 & model_est{2}(:,5)==nperm ,1);
                tstat_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)>0 & model_est{2}(:,5)==nperm ,3);
            else
                sig_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)<0 & model_est{2}(:,5)==nperm ,1);
                tstat_elec=model_est{2}(model_est{2}(:,4)<clus_alpha & model_est{2}(:,3)<0 & model_est{2}(:,5)==nperm ,3);
            end
            if ~isempty(sig_elec)
                perm_clus{nsign,nperm}{1,1}={neighbours(sig_elec(1)).label};
                perm_clus{nsign,nperm}{1,3}=tstat_elec(1);
                perm_clus{nsign,nperm}{1,2}=neighbours(sig_elec(1)).neighblabel';
                for nEl=2:length(sig_elec)
                    clus_found=0;
                    for nclus=1:size(perm_clus{nsign,nperm},1)
                        if sum(ismember(perm_clus{nsign,nperm}{nclus,1},neighbours(sig_elec(nEl)).neighblabel'))~=0
                            perm_clus{nsign,nperm}{nclus,1}=unique([perm_clus{nsign,nperm}{nclus,1} ; {neighbours(sig_elec(nEl)).label}]);
                            perm_clus{nsign,nperm}{nclus,3}=[perm_clus{nsign,nperm}{nclus,3} ; tstat_elec(nEl)];
                            perm_clus{nsign,nperm}{nclus,2}=unique([perm_clus{nsign,nperm}{nclus,2} ; neighbours(sig_elec(nEl)).neighblabel']);
                            clus_found=1;
                            break;
                        end
                    end
                    if clus_found==0
                        num_clus=size(perm_clus{nsign,nperm},1);
                        perm_clus{nsign,nperm}{num_clus+1,1}= {neighbours(sig_elec(nEl)).label};
                        perm_clus{nsign,nperm}{num_clus+1,3}=tstat_elec(nEl);
                        perm_clus{nsign,nperm}{num_clus+1,2}=neighbours(sig_elec(nEl)).neighblabel';
                    end
                end
                if ~isempty(limSize)
                    clus_discard=[];
                    for m=1:size(perm_clus{nsign,nperm},1)
                        if size(perm_clus{nsign,nperm}{m,1},1)<=limSize
                            clus_discard=[clus_discard m];
                        end
                    end
                    perm_clus{nsign,nperm}(clus_discard,:)=[];
                end
                clus_merged=[];
                for nclus=1:size(perm_clus{nsign,nperm},1)
                    if ismember(nclus,clus_merged)
                        continue;
                    end
                    for nclus2=nclus+1:size(perm_clus{nsign,nperm},1)
                        if sum(ismember(perm_clus{nsign,nperm}{nclus2,1},perm_clus{nsign,nperm}{nclus,2}))~=0
                            perm_clus{nsign,nperm}{nclus,3}=[perm_clus{nsign,nperm}{nclus,3} ; perm_clus{nsign,nperm}{nclus2,3}(~ismember(perm_clus{nsign,nperm}{nclus2,1},perm_clus{nsign,nperm}{nclus,1}))];
                            perm_clus{nsign,nperm}{nclus,1}=unique([perm_clus{nsign,nperm}{nclus,1} ; perm_clus{nsign,nperm}{nclus2,1}]);
                            perm_clus{nsign,nperm}{nclus,2}=unique([perm_clus{nsign,nperm}{nclus,2} ; perm_clus{nsign,nperm}{nclus2,2}]);
                            clus_merged=[clus_merged nclus2];
                        end% Model_est = Elec x ? x t value x p value x group?

                    end
                end
                perm_clus{nsign,nperm}(clus_merged,:)=[];                
            end
        end
    end
    
    perm_tval=cell(1,2);
    for nsign=1:2
        for nperm=1:totperm
            temp_tval=[];
            for nclus=1:size(perm_clus{nsign,nperm},1)
                temp_tval=[temp_tval sum(perm_clus{nsign,nperm}{nclus,3})];
            end
            if nsign==1
                if isempty(temp_tval)
                    perm_tval{nsign}=[perm_tval{nsign} 0];
                else
                    perm_tval{nsign}=[perm_tval{nsign} max(temp_tval)];
                end
            else
                if isempty(temp_tval)
                    perm_tval{nsign}=[perm_tval{nsign} 0];
                else
                    perm_tval{nsign}=[perm_tval{nsign} min(temp_tval)];
                end
            end
        end
    end
    
    all_clus=[]; nc=0;
    for nsign=1:2
        for nclus=1:size(real_clus{nsign},1)
            this_tval=sum(real_clus{nsign}{nclus,3});
            if nsign==1
                if mean(this_tval<perm_tval{nsign})<montecarlo_alpha
                    nc=nc+1;
                    all_clus{nc}={'pos' real_clus{nsign}{nclus,1} this_tval sum(this_tval<perm_tval{nsign})/totperm};
                end
            else
                if mean(this_tval>perm_tval{nsign})<montecarlo_alpha
                    nc=nc+1;
                    all_clus{nc}={'neg' real_clus{nsign}{nclus,1} this_tval sum(this_tval>perm_tval{nsign})/totperm};
                end
            end
        end
    end
end
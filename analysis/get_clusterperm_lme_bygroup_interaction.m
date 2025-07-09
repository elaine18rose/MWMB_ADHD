function [all_clus]=get_clusterperm_lme_bygroup_interaction(model_est,clus_alpha,montecarlo_alpha,totperm,neighbours,limSize)

if nargin<6
    limSize=[]; % ???
end
for nEffect=1:2 % Main effect and interaction
    real_clus=cell(1,2);
    for nsign=1
        sig_elec=model_est{1}(model_est{1}(:,1+2*nEffect)<clus_alpha,1);
        Fstat_elec=model_est{1}(model_est{1}(:,1+2*nEffect)<clus_alpha,2*nEffect);
        if ~isempty(sig_elec)
            clusters = find_clusters(sig_elec, neighbours);
            for nclus=1:length(clusters)
                real_clus{nsign}{nclus,1}=clusters{nclus}; % Get the names of its neighbours
                real_clus{nsign}{nclus,3}=sum(Fstat_elec(ismember({neighbours(sig_elec).label},clusters{nclus}))); % And its t value
                real_clus{nsign}{nclus,2}=[]; % And its neighbours? (how is this different?)
            end
        end
    end

    perm_clus=cell(2,totperm);
    for nperm=1:totperm
        for nsign=1
            sig_elec=model_est{2}(model_est{2}(:,1+2*nEffect)<clus_alpha & model_est{2}(:,6)==nperm,1);
            Fstat_elec=model_est{2}(model_est{2}(:,1+2*nEffect)<clus_alpha & model_est{2}(:,6)==nperm,2*nEffect);

            if ~isempty(sig_elec)
                clusters = find_clusters(sig_elec, neighbours);
                for nclus=1:length(clusters)
                    perm_clus{nsign,nperm}{nclus,1}=clusters{nclus}; % Get the names of its neighbours
                    perm_clus{nsign,nperm}{nclus,3}=sum(Fstat_elec(ismember({neighbours(sig_elec).label},clusters{nclus}))); % And its t value
                    perm_clus{nsign,nperm}{nclus,2}=[]; % And its neighbours? (how is this different?)
                end
            end
        end
    end

    perm_Fval=cell(1,2);
    for nsign=1 %:2
        for nperm=1:totperm
            temp_Fval=[];
            for nclus=1:size(perm_clus{nsign,nperm},1)
                temp_Fval=[temp_Fval sum(perm_clus{nsign,nperm}{nclus,3})];
            end
            if isempty(temp_Fval)
                perm_Fval{nsign}=[perm_Fval{nsign} 0];
            else
                perm_Fval{nsign}=[perm_Fval{nsign} max(temp_Fval)];
            end
        end
    end

    all_clus{nEffect}=[]; nc=0;
    for nsign=1 %:2
        for nclus=1:size(real_clus{nsign},1)
            this_Fval=sum(real_clus{nsign}{nclus,3});
            if mean(this_Fval>perm_Fval{nsign})>1-montecarlo_alpha
                nc=nc+1;
                all_clus{nEffect}{nc}={'pos' real_clus{nsign}{nclus,1} this_Fval sum(this_Fval<perm_Fval{nsign})/totperm};
            end
        end
    end
end

load('/Users/thandrillon/Downloads/all_erp.mat')
totSubj=size(all_erp{1},1);
vec_numSubj=floor([totSubj totSubj/2 totSubj/4 totSubj/8 totSubj/16 totSubj/32 totSubj/64]);
montecarloalpha=0.05;
clusteralpha=0.1;
npermutation=100;
States={'W','1','2','3','R'};
figure;
for k=1:5
    data=[];
    data{1}(1,:,:)=squeeze(all_erp{k}(~isnan(all_erp{k}(:,1,1)),1,:))';

    [realpos realneg]=get_cluster_permutation(data,montecarloalpha,clusteralpha,npermutation,1:size(data{1},2),1);
    posclusterids=find(realpos{1}.pmonte<montecarloalpha);
    negclusterids=find(realneg{1}.pmonte<montecarloalpha);
    lineclusters{1,k}=double(ismember(realpos{1}.clusters,posclusterids))-1*double(ismember(realneg{1}.clusters,negclusterids));
    subplot(length(vec_numSubj),5,k);
    imagesc(lineclusters{1,k});
    title(sprintf('%s - N=%g',States{k},vec_numSubj(1)))
    caxis([-1 1])
end

for m=2:length(vec_numSubj)
    subsample=100;
    for k=1:5
        for j=1:subsample
            data=[];
            data{1}(1,:,:)=squeeze(all_erp{k}(~isnan(all_erp{k}(:,1,1)),1,:))';
            subids=randperm(size(data{1},3));
            subids=subids(1:vec_numSubj(m));
            data{1}= data{1}(1,:,subids);

            [realpos realneg]=get_cluster_permutation(data,montecarloalpha,clusteralpha,npermutation,1:size(data{1},2),1);
            posclusterids=find(realpos{1}.pmonte<montecarloalpha);
            negclusterids=find(realneg{1}.pmonte<montecarloalpha);
            lineclusters{m,k}(j,:)=double(ismember(realpos{1}.clusters,posclusterids))-1*double(ismember(realneg{1}.clusters,negclusterids));
        end
        subplot(length(vec_numSubj),5,k+(m-1)*5);
        imagesc(lineclusters{m,k});
        title(sprintf('%s - N=%g',States{k},vec_numSubj(m)))
        caxis([-1 1])
    end

end
% format_fig;

%%
figure;
for k=1:5
    subplot(1,5,k)
    for m=1:length(vec_numSubj)

        hold on;
        plot(mean(lineclusters{m,k},1))
    end
    line(xlim,[1 1]*0.8,'Color','k','LineWidth',2)
    line(xlim,[1 1]*-0.8,'Color','k','LineWidth',2)
    format_fig;
    if k==5
        legend({'729','364','182','91','45','22','11'})
    end
    title(States{k})
end
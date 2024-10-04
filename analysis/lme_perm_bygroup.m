function [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_bygroup(table,predictor,formula,totperm,all_pred_perm)
% table=GO_table;
%
% formula='SW~1+Group+(1|SubID)';
% permvars={'SubID'};
% totperm=100;
if nargin<5
    all_pred_perm=[];
end
out_pred_perm=[];
% run real model
eval(sprintf('table.pred=table.%s;',predictor));
model= fitlme(table,formula);
real_out=[double(model.Coefficients(find_trials([model.CoefficientNames],'pred_'),2)) double(model.Coefficients(find_trials([model.CoefficientNames],'pred_'),4)) double(model.Coefficients(find_trials([model.CoefficientNames],'pred_'),6))];
cont_out=[model.CoefficientNames(find_trials([model.CoefficientNames],'pred_'))];

uniqueIDs=unique(table.SubID);
for nSub=1:length(uniqueIDs)
    uniqueGroups(nSub)=unique(table.(predictor)(table.SubID==uniqueIDs(nSub)));
end
perm_out=[];
cont_perm_out=[];
fprintf('%4.0f/%4.0f\n',0,totperm)
for np=1:totperm
    if size(all_pred_perm,1)~=totperm
        group_perm_idx=1:length(uniqueGroups);
        group_perm_idx=group_perm_idx(randperm(length(group_perm_idx)));
        out_pred_perm(np,:)=group_perm_idx;
    else
        group_perm_idx=all_pred_perm(np,:)';
        out_pred_perm(np,:)=group_perm_idx;
    end
    uniqueGroups_perm=uniqueGroups(group_perm_idx);
    table_perm=table;
    for nS=1:length(uniqueIDs)
        table_perm.pred(table.SubID==uniqueIDs(nS))=uniqueGroups_perm(nS);
    end
    model_perm= fitlme(table_perm,formula);
    ContNames=unique(table_perm.pred);
    perm_out=[perm_out ; double(model_perm.Coefficients(find_trials([model_perm.CoefficientNames],'pred_'),2)) double(model_perm.Coefficients(find_trials([model_perm.CoefficientNames],'pred_'),4)) double(model_perm.Coefficients(find_trials([model_perm.CoefficientNames],'pred_'),6)) np];
    cont_perm_out=[cont_perm_out ; model_perm.CoefficientNames(find_trials([model_perm.CoefficientNames],'pred_'))];
      fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
end
fprintf('\n');

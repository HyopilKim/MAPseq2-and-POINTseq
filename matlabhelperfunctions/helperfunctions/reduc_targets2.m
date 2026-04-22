function [B_reduc,Bseq_reduc,spikes_reduc]=reduc_targets(barcodematrix,spikes)

%reduced Bnorm
Bnorm_input=[];
for i=1:size(BnormCell,2)
    Bnorm_input=[Bnorm_input;BnormCell{i}];
end
Bnorm_reduc=zeros(size(Bnorm_input,1),size(reducdim,2)+1);
for i=1:size(reducdim,2)
    Bnorm_reduc(:,i)=sum(Bnorm_input(:,reducdim{i}),2);
end

Bnorm_reduc(:,size(reducdim,2)+1)=Bnorm_input(:,end);
Bnorm_reduc=Bnorm_reduc(max(Bnorm_reduc(:,1:size(reducdim,2)),[],2)>0,:);
maxBnorm_reduc=Bnorm_reduc(:,1:size(reducdim,2))./repmat(max(Bnorm_reduc(:,1:size(reducdim,2)),[],2),1,size(Bnorm_reduc(:,1:size(reducdim,2)),2));
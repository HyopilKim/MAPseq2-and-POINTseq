function [barcodematrix89,spikes89]=reduc_targets(barcodematrix,spikes,old,new,sorting95to89)

%reduced Bnorm
Bnorm_input=[];
for i=1:size(BnormCell,2)
    Bnorm_input=[Bnorm_input;BnormCell{i}];
end
Bnorm_reduc=zeros(size(Bnorm_input,1),size(reducdim,2)+1);
for i=1:size(reducdim,2)
    Bnorm_reduc(:,i)=sum(Bnorm_input(:,reducdim{i}),2);
end

Bnorm_reduc(:,size(reducdim,2)+1)=max(Bnorm_input(:,sourcesite),[],2);
Bnorm_reduc=Bnorm_reduc(max(Bnorm_reduc(:,1:size(reducdim,2)),[],2)>0,:);
maxBnorm_reduc=Bnorm_reduc(:,1:size(reducdim,2))./repmat(max(Bnorm_reduc(:,1:size(reducdim,2)),[],2),1,size(Bnorm_reduc(:,1:size(reducdim,2)),2));

barcodematrix89=barcodematrix(:,sorting95to89);
spikes89=spikes(sorting95to89);

for i=1:size(old,2)
    for j=1:length(old{i})
    barcodematrix89_dat8(:,new{i})=barcodematrix_dat8(:,old{j})

    
    barcodematrix89_dat8(:,6)=barcodematrix_dat8(:,6)+barcodematrix(:,8);
barcodematrix89_dat8(:,81)=barcodematrix(:,85)+barcodematrix(:,91);
barcodematrix89(:,82)=barcodematrix(:,86)+barcodematrix(:,92);


spikes89=combine_spikes(old,new,spikes,spikes89);


spikes89=combine_spikes([6 8],6,spikes,spikes89);
spikes89=combine_spikes([85 91],81,spikes,spikes89);
spikes89=combine_spikes([86 92],82,spikes,spikes89);
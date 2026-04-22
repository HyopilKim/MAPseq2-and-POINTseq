function [barcodematrix_reduc,spikes_reduc]=reduc_targets(barcodematrix,spikes,reducdim,sourcesite)

%reduced barcodematrix and spikes
barcodematrix_reduc=zeros(size(barcodematrix,1),size(reducdim,2));
spikes_reduc=struct();
for i=1:size(reducdim,2)
    barcodematrix_reduc(:,i)=sum(barcodematrix(:,reducdim{i}),2);
    spikes_reduc(i).counts2u=[];
    for j=1:length(reducdim{i})
        spikes_reduc(i).counts2u=[spikes_reduc(i).counts2u;spikes(reducdim{i}(j)).counts2u];
    end
end
barcodematrix_reduc(:,size(reducdim,2)+1)=barcodematrix(:,sourcesite)
spikes_reduc(size(reducdim,2)+1).counts2u=spikes(sourcesite).counts2u;

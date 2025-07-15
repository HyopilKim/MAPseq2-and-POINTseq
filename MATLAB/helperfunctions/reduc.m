function [barcodematrix_reduc,refbarcodes_reduc,spikes_reduc]=reduc(barcodematrix,refbarcodes,spikes,sorting,old,new,filter1,filter2)


barcodematrix_reduc=barcodematrix(:,sorting);




if exist('old','var') && exist('new','var')
    for i=1:size(old,2)
        sum=zeros(size(barcodematrix,1),1);
        dim=old{i};
        for j=1:length(old{i})
            sum=sum+barcodematrix(:,dim(j));
        end
     barcodematrix_reduc(:,new(i))=sum;
    end
    spikes_reduc=combine_spikes(old,new,spikes,spikes(sorting));
else
    spikes_reduc=spikes(sorting);
end

if exist('filter1','var') && exist('filter2','var')
    barcodematrix_reduc=barcodematrix_reduc(filter1,:);
    barcodematrix_reduc=barcodematrix_reduc(filter2,:);
end

nonzero=max(barcodematrix_reduc,[],2)>0;
barcodematrix_reduc=barcodematrix_reduc(nonzero,:);
refbarcodes_reduc=refbarcodes(nonzero,:);



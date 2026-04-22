function [B,Bnorm]=normBCmat1(barcodematrix,spikes,sourcethresh,projthresh,sourcesite,projsite,cellbodythresh)


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end

%filter barcodes raw UMI by thresholds
if exist('sourcethresh','var') && exist('projthresh','var') && exist('sourcesite','var') && exist('projsite','var')
   B=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
else
   B=barcodematrix;
end

%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);

%filter by cellbodythresh
if exist('cellbodythresh','var')
    B=B(max(B_tar,[],2)<cellbodythresh,:);
    Bnorm=Bnorm(max(B_tar,[],2)<cellbodythresh,:);
end

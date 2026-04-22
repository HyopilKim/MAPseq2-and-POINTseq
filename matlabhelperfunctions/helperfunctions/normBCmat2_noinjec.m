function [B,Bnorm,Bseq]=normBCmat2_noinjec(barcodematrix,refbarcodes,spikes,projthresh)


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end


%filter barcodes raw UMI by thresholds
B=barcodematrix(max(barcodematrix,[],2)>projthresh,:);
Bseq=refbarcodes(max(barcodematrix,[],2)>projthresh,:);

%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);


save('filtBCmat.mat','B','Bnorm','Bseq');

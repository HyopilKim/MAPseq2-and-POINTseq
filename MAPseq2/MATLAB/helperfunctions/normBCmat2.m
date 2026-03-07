function normBCmat2(barcodematrix,refbarcodes,spikes,sourcethresh,projthresh,sourcesite,projsite,sorting,cellbodythresh)

%sourcethresh and projthresh are thresholds for source and target regions,respectively
%sorting is reordering target regions
% cellbodythresh is a limit for target region UMI to exclude potential propagated barcodes


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end


%filter barcodes raw UMI by thresholds
B=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
Bseq=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);

%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);

%targets only
B_tar=B(:,projsite);
Bnorm_tar=Bnorm(:,projsite);

if exist('sorting','var')
    B_tar=B(:,sorting);
    Bnorm_tar=Bnorm(:,sorting);
end


%filter by cellbodythresh
if exist('cellbodythresh','var')
    B=B(max(B_tar,[],2)<cellbodythresh,:);
    B_tar=B_tar(max(B_tar,[],2)<cellbodythresh,:);
    Bnorm=Bnorm(max(B_tar,[],2)<cellbodythresh,:);
    Bnorm_tar=Bnorm_tar(max(B_tar,[],2)<cellbodythresh,:);
    Bseq=Bseq(max(B_tar,[],2)<cellbodythresh,:);
end

%normalization by total projection or max projection
totalBnorm_tar=Bnorm_tar./sum(Bnorm_tar,2);
maxBnorm_tar=Bnorm_tar./repmat(max(Bnorm_tar,[],2),1,size(Bnorm_tar,2));

save('filtBCmat.mat','B','Bnorm','Bseq','B_tar','Bnorm_tar','totalBnorm_tar','maxBnorm_tar');

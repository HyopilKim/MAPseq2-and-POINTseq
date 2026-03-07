function normBCmat3(barcodematrix,refbarcodes,spikes,sourcethresh,projthresh,sourcesite,projsite,ratiothreshold,sorting,sourcelimit)


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end


%filter barcodes raw UMI by thresholds
B=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
Bseq=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);

%filter by cellbodythresh and sourcelimit
if exist('sourcelimit','var') 
    B=B(max(B(:,sourcesite),[],2)<sourcelimit,:);
    Bseq=Bseq(max(B(:,sourcesite),[],2)<sourcelimit,:);
end


%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);

%targets only
B_tar=B(:,projsite);
Bnorm_tar=Bnorm(:,projsite);


if exist('sorting','var')
    B_tar=B(:,sorting);
    Bnorm_tar=Bnorm(:,sorting);
end

%normalization by total projection or max projection
maxBnorm_tar=Bnorm_tar./repmat(max(Bnorm_tar,[],2),1,size(Bnorm_tar,2));
totalBnorm_tar=Bnorm_tar./sum(Bnorm_tar,2);

%separate by injsite (in Cell variables, last column indicates the sourcesite where 1=SNc, 2=VTA) 
if exist('ratiothreshold','var')
    Bnorm_source=Bnorm(:,sourcesite);
    [sourcemax,index]=max(Bnorm_source,[],2);
    ratio=sort(Bnorm_source./sourcemax,2,'descend');
    B=B(ratio(:,2)<ratiothreshold,:);
    Bseq=Bseq(ratio(:,2)<ratiothreshold,:);
    Bnorm=Bnorm(ratio(:,2)<ratiothreshold,:);
    B_tar=B_tar(ratio(:,2)<ratiothreshold,:);
    Bnorm_tar=Bnorm_tar(ratio(:,2)<ratiothreshold,:);
    maxBnorm_tar=maxBnorm_tar(ratio(:,2)<ratiothreshold,:);
    totalBnorm_tar=totalBnorm_tar(ratio(:,2)<ratiothreshold,:);
    index=index(ratio(:,2)<ratiothreshold);
    filter2=(ratio(:,2)<ratiothreshold);
else
    index=[];
    filter2=[];
end

filter1=(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh);


save('filtBC.mat','B_tar','Bnorm_tar','B','Bnorm','Bseq','totalBnorm_tar','maxBnorm_tar','index','filter1','filter2');
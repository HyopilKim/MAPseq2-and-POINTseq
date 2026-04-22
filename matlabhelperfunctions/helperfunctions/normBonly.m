function [Bnorm,BCell,BnormCell,BseqCell,Bnorm_tarCell,maxBnorm_tarCell,totalBnorm_tarCell]=normBCmat_sources(barcodematrix,refbarcodes,spikes,sourcethresh,projthresh,sourcesite,projsite,sorting,ratiothreshold)


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
Bnorm_tar=Bnorm_tar(:,sorting);
totalBnorm_tar=Bnorm_tar./sum(Bnorm_tar,2)
maxBnorm_tar=Bnorm_tar./repmat(max(Bnorm_tar,[],2),1,size(Bnorm_tar,2))

%separate by injsite
Bnorm_source=Bnorm(:,sourcesite);
[sourcemax,index]=max(Bnorm_source,[],2);
ratio=sort(Bnorm_source./sourcemax,2,'descend');

BseqCell={};
for i=1:length(sourcesite)
    Bseq=double(Bseq)
    BseqCell{i}=Bseq(ratio(:,2)<ratiothreshold & index==i,:);
    BseqCell{i}=[BseqCell{i},repmat(i,size(BseqCell{i},1),1)];
end

BCell={};
for i=1:length(sourcesite)
    BCell{i}=B(ratio(:,2)<ratiothreshold & index==i,:);
    BCell{i}=[BCell{i},repmat(i,size(BCell{i},1),1)];
end

BnormCell={};
for i=1:length(sourcesite)
    BnormCell{i}=Bnorm(ratio(:,2)<ratiothreshold & index==i,:);
    BnormCell{i}=[BnormCell{i},repmat(i,size(BnormCell{i},1),1)];
end

Bnorm_tarCell={};
for i=1:length(sourcesite)
    Bnorm_tarCell{i}=Bnorm_tar(ratio(:,2)<ratiothreshold & index==i,:);
    Bnorm_tarCell{i}=[Bnorm_tarCell{i},repmat(i,size(Bnorm_tarCell{i},1),1)];
end

maxBnorm_tarCell={};
for i=1:length(sourcesite)
    maxBnorm_tarCell{i}=maxBnorm_tar(ratio(:,2)<ratiothreshold & index==i,:);
    maxBnorm_tarCell{i}=[maxBnorm_tarCell{i},repmat(i,size(maxBnorm_tarCell{i},1),1)];
end

totalBnorm_tarCell={};
for i=1:length(sourcesite)
    totalBnorm_tarCell{i}=totalBnorm_tar(ratio(:,2)<ratiothreshold & index==i,:);
    totalBnorm_tarCell{i}=[totalBnorm_tarCell{i},repmat(i,size(totalBnorm_tarCell{i},1),1)];
end

save('filtBCmat.mat','BCell','BseqCell','BnormCell','Bnorm_tarCell','maxBnorm_tarCell','totalBnorm_tarCell');

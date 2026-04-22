function Bnorm=normBCmat0(barcodematrix,spikes)


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end


B=barcodematrix;


%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);


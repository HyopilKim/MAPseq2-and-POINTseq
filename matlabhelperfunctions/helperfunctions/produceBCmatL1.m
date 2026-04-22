function produceBCmatL1(bcn,prefix,brain)
%process MAPseq data to produce barcode matrix

load rawdata

%bcn=1:length(bcn);

%filter out SSI with no barcodes
filesize=[];
for i=1:length(bcn)
    file=dir(['barcodesplitter/thresholds/',prefix,int2str(bcn(i)),'_countsL1.txt']);
    filesize(i)=file.bytes;
end

bcnfilt=bcn(filesize>0);
data=struct();
for i=1:length(bcnfilt)
data(i).counts=dlmread(['barcodesplitter/thresholds/',prefix,int2str(bcnfilt(i)),'_countsL1.txt']);
data(i).reads=int8(char(textread(['barcodesplitter/thresholds/',prefix,int2str(bcnfilt(i)),'_seqL1.txt'],'%s')));
end
 
 save('dataL1.mat','data','-v7.3');
 
 
 
 
%% Finish error correction by reading in bowtie alignments and finding connected graph components
 
%load alignments
positions1=[];
for i=1:length(bcnfilt)
positions1(i).x=dlmread(['barcodesplitter/thresholds/','bowtie',int2str(bcnfilt(i)),'_2u_1L1.txt']);
positions1(i).y=dlmread(['barcodesplitter/thresholds/','bowtie',int2str(bcnfilt(i)),'_2u_3L1.txt']);
clustermatrix1(i).C=sparse(positions1(i).x,positions1(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
end
 %save('clustermatrix1L1.mat','clustermatrix1','-v7.3');
% 
 %load clustermatrix1
 
%find connected components
graph=[];
for i=1:length(bcnfilt)
i
    [graph(i).S,graph(i).G]=graphconncomp(clustermatrix1(i).C,'Directed','false'); %find the connected graph components
end
% save('graph1L1.mat','graph');
% 
% load graph1
 
%collapse barcodes to most abundant member of the connected graph component
 
for i=1:length(bcnfilt)
x=1:graph(i).S;
[tf,loc]=ismember(x,graph(i).G,'R2012a');
collapsedreads=data(i).reads(loc,:);
collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
[corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
corrected(i).reads2u=collapsedreads(ix,:);
end
 
 
%% remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
for i=1:length(bcnfilt)
    a=findhomopolymers(corrected(i).reads2u,minrunlength);
    corrected(i).freads=corrected(i).reads2u(~a,:);
    corrected(i).fcounts=corrected(i).counts2u(~a,:);
end
 
for i=1:length(bcnfilt)
data(i).BCseqff=corrected(i).freads;
data(i).BCcountsff=corrected(i).fcounts;
end
save('dataL1.mat','data','-v7.3');
 
 
 %% load data.m and match the barcodes from the different areas against each other to produce a barcode matrix
% this is code modified to in situ MAPseq, so it will run without an
% injection site.


%% find overlap of barcodes
%% set up reference set of barcodes
%load data
load dataL1.mat

%collect all suqeunces detected in the target sites
refbarcodes_tmp=[];
for i=1:length(data)
refbarcodes_tmp=[refbarcodes_tmp;data(i).BCseqff];  
end
refbarcodes=unique(refbarcodes_tmp,'rows'); %unique barocodes to get reference set

%% construct barcodematrix by matching barcodes in target sites to reference barcodes

barcodematrix=zeros(size(refbarcodes,1),length(data));%initiate the barcode matrix

for i=1:length(data)
    %pull out reads and counts into new variables for ease of use
    BCreads=data(i).BCseqff;
    BCmolcounts=data(i).BCcountsff;
    [ind,loc]=ismember(BCreads,refbarcodes,'rows');
    barcodematrix(loc(loc~=0),i)=BCmolcounts(ind);   
end


%set thresholds of how many molecule counts a barcode has to have to be
%considered real. by default we leave this at 0
%threshold=0; %lowest molecule count considered trustworthy
%barcodematrixthresholded=barcodematrix(sum(barcodematrix>threshold,2)~=0,:);
%refbarcodesthresholded=refbarcodes(sum(barcodematrix>threshold,2)~=0,:);
                           barcodematrixL1 = barcodematrix;
                           refbarcodesL1 = refbarcodes;
                           tmp = zeros(size(barcodematrix,1),length(bcn));
                           tmp(:,bcnfilt) = barcodematrixL1;
                           barcodematrixL1 = tmp;

save(['barcodematrixL1',prefix,'_',brain,'.mat'],'barcodematrixL1','refbarcodesL1','-v7.3')

						   

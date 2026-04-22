function globalcollapse(prefix,bcn,spikeseq,brain)
%to provide globally collapsed refbarcodes for bowtie
filesize=[];
for i=1:length(bcn)
    file=dir([prefix,'_',int2str(bcn(i)),'_counts.txt']);
    filesize(i)=file.bytes;
end

bcnfilt=bcn(filesize>0);

data=struct();
data.counts=dlmread([prefix,'_global_counts.txt']);
data.reads=int8(char(textread([prefix,'_global_seq.txt'],'%s')));
for i=1:size(data.reads,1)
    isspikes(i)=all(data.reads(i,25:32)==int8(spikeseq));
end

%save(['data_',brain,'.mat','data','-v7.3');
  
 
%% Finish error correction by reading in bowtie alignments and finding connected graph components
 
%load alignments
positions1=[];
positions1.x=dlmread('bowtieglobal_2u_1.txt');
positions1.y=dlmread('bowtieglobal_2u_3.txt');
clustermatrix1.C=sparse(positions1.x,positions1.y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries

 %save('clustermatrix1.mat','clustermatrix1','-v7.3');
% 
 %load clustermatrix1
 
%find connected components
graph=[];
[graph.S,graph.G]=graphconncomp(clustermatrix1.C,'Directed',false); %find the connected graph components
% save('graph1.mat','graph');
% 
% load graph1
 
%collapse barcodes to most abundant member of the connected graph component
G=sort(unique(graph.G));
collapsedreadsall=data.reads(graph.G,:);
G_spikes=sort(unique(graph.G(isspikes)));
G_nospikes=sort(setdiff(G, G_spikes));
[~,loc_spikes]=ismember(G_spikes,graph.G,'R2012a');
[~,loc_nospikes]=ismember(G_nospikes,graph.G,'R2012a');
collapsedreads_spikes=data.reads(loc_spikes,:);
collapsedreads_nospikes=data.reads(loc_nospikes,:);


%% remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
a_spikes=findhomopolymers(collapsedreads_spikes,minrunlength);
corrected_spikes=collapsedreads_spikes(~a_spikes,:);
a_nospikes=findhomopolymers(collapsedreads_nospikes,minrunlength);
corrected_nospikes=collapsedreads_nospikes(~a_nospikes,:);

seqs=struct();
for i=1:size(bcnfilt,2)
    seqs(i).seq=int8(char(textread(['DAT14_',int2str(bcnfilt(i)),'_seq.txt'],'%s')));
    seqs(i).count=textread(['DAT14_',int2str(bcnfilt(i)),'_counts.txt']);
end

spikes=struct();
for i = 1:numel(seqs)
    eachseq=int8(char(seqs(i).seq));
    [~,loc_each]=ismember(eachseq,data.reads,'rows');
    [eachcollapsedseq,~,grp]=unique(collapsedreadsall(loc_each,:),'rows','stable');
    eachcollapsedcount=accumarray(grp, seqs(i).count);
    seqs(i).collapsedseq = eachcollapsedseq;
    seqs(i).collapsedcount  = eachcollapsedcount;

    tf_eachspikes=ismember(seqs(i).collapsedseq,corrected_spikes,'rows');
    spikes(i).reads2u=seqs(i).collapsedseq(tf_eachspikes,:);
    spikes(i).counts2u=seqs(i).collapsedcount(tf_eachspikes,:);
    
    tf_eachnospikes=ismember(seqs(i).collapsedseq,corrected_nospikes,'rows');
    seqs(i).seq_nospikes=seqs(i).collapsedseq(tf_eachnospikes,:);
    seqs(i).count_nospikes=seqs(i).collapsedcount(tf_eachnospikes,:);
end
    


%collect all suquences detected in the target sites
refbarcodes_tmp=[];
for i=1:length(seqs)
refbarcodes_tmp=[refbarcodes_tmp;seqs(i).seq_nospikes];  
end
refbarcodes=unique(refbarcodes_tmp,'rows'); %unique barocodes to get reference set

%% construct barcodematrix by matching barcodes in target sites to reference barcodes

barcodematrix=zeros(size(refbarcodes,1),length(seqs));%initiate the barcode matrix

for i=1:length(data)
    %pull out reads and counts into new variables for ease of use
    BCreads=seqs(i).seq_nospikes;
    BCmolcounts=seqs(i).count_nospikes;
    [ind,lc]=ismember(BCreads,refbarcodes,'rows');
    barcodematrix(lc(lc~=0),i)=BCmolcounts(ind);   
end

%set thresholds of how many molecule counts a barcode has to have to be
%considered real. by default we leave this at 0
threshold=0; %lowest molecule count considered trustworthy
barcodematrix=barcodematrix(sum(barcodematrix>threshold,2)~=0,:);
refbarcodes=refbarcodes(sum(barcodematrix>threshold,2)~=0,:);


save(['barcodematrix',prefix,'_',brain,'.mat'],'barcodematrix','refbarcodes','-v7.3')
save(['spikes',prefix,'_',brain,'.mat'],'spikes')
save(['seqs',prefix,'_',brain,'.mat'],'seqs')

						   
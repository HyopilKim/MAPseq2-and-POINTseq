function globalcollapse(prefix,spikeseq)
%to provide globally collapsed refbarcodes for bowtie


data=struct();
data.counts=dlmread(['thresholds/',prefix,'_global_counts.txt']);
data.reads=int8(char(textread(['thresholds/',prefix,'_global_seq.txt'],'%s')));
for i=1:size(data.reads,1)
    isspikes(i)=all(data.reads(i,25:32)==int8(spikeseq));
end

%filter out SSI with no barcodes
filesize=[];
for i=1:length(bcn)
    file=dir(['thresholds/',prefix,'_',int2str(bcn(i)),'_counts.txt']);
    filesize(i)=file.bytes;
end

bcnfilt=bcn(filesize>0);
data=struct();
for i=1:length(bcnfilt)
data(i).counts=dlmread(['thresholds/',prefix,'_',int2str(bcnfilt(i)),'_counts.txt']);
data(i).reads=int8(char(textread(['thresholds/',prefix,'_',int2str(bcnfilt(i)),'_seq.txt'],'%s')));
end




%load alignments
positions1=[];
positions1.x=dlmread(['thresholds/','bowtieglobal_2u_1.txt']);
positions1.y=dlmread(['thresholds/','bowtieglobal_2u_3.txt']);
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
G=unique(graph.G);
G_spikes=sort(unique(graph.G(isspikes)));
[~,loc_spikes]=ismember(G_spikes,graph.G,'R2012a');
collapsedreads_spikes=data.reads(loc_spikes,:);
collapsedreads_nospikes=data.reads(~loc_spikes,:);
  
%% remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
a_spikes=findhomopolymers(collapsedreads_spikes,minrunlength);
corrected_spikes=collapsedreads_spikes(~a_spikes,:);
a_nospikes=findhomopolymers(collapsedreads_nospikes,minrunlength);
corrected_nospikes=collapsedreads_nospikes(~a_nospikes,:);

refspikes=string(char(corrected_spikes));
writelines(refspikes, "refspikes.txt");
refnospikes=string(char(corrected_nospikes));
writelines(refnospikes, "refnospikes.txt");

save(['globalcollapsed',prefix,'.mat'],'refspikes','refnospikes','-v7.3')

						   
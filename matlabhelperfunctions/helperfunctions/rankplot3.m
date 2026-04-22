function rankplot3(bcn,prefix,donotsave)
%plot rank plot of counts
%% read in sequences and counts of data split into the different libraries

%bcn=212:229;
for i=1:length(bcn)
data(i).rawcounts=dlmread([prefix,'_BC',int2str(bcn(i)),'_BU_counts.txt']);
end
save('rawdata.mat','data','-v7.3');

%use this to look at the sequence rank plot of every libary and choose a threshold for preprocessing.sh 
load rawdata
if ~exist('donotsave','var')
    figure;for i=1:length(bcn);semilogy(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');
        ylabel('Read count');findfigs;exportgraphics(gcf, sprintf('figure_%d.png', bcn(i)));end
else
    figure;for i=1:length(bcn);semilogy(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');
        ylabel('Read count');findfigs;pause;end
end


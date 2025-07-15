function rankplot2(bcn,prefix,shape,xlimit)
%unless needed, leave shape as 0 to have semilogy
%plot rank plot of counts
%% read in sequences and counts of data split into the different libraries

%bcn=212:229;
for i=1:length(bcn)
data(i).rawcounts=dlmread([prefix,'_BC',int2str(bcn(i)),'_BU_counts.txt']);
end
save('rawdata.mat','data','-v7.3');

%use this to look at the sequence rank plot of every libary and choose a threshold for preprocessing.sh 
load rawdata
if shape==0 && exist('xlimit','var')
    figure;for i=1:length(bcn);semilogy(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');xlim([0 xlimit]);
        ylabel('Read count');findfigs;exportgraphics(gcf, sprintf('figure_%d.png', i));end
elseif shape==1 && exist('xlimit','var')
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');xlim([0 xlimit]);
        ylabel('Read count');findfigs;exportgraphics(gcf, sprintf('figure_%d.png', i));end
elseif shape==0 && ~exist('xlimit','var')
    figure;for i=1:length(bcn);semilogy(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');
        ylabel('Read count');findfigs;exportgraphics(gcf, sprintf('figure_%d.png', i));end
elseif shape==1 && ~exist('xlimit','var')
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(bcn(i)));xlabel('Sequence rank');
        ylabel('Read count');findfigs;exportgraphics(gcf, sprintf('figure_%d.png', i));end
end


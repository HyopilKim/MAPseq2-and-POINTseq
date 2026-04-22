function [input1,input2]=sepclustergram(matrix,targets,sources)

input=matrix(:,1:(end-1));
idx=matrix(:,end);

cg=clustergram(input,'Cluster','column','Symmetric','false','Colormap','parula','ColumnLabels',targets);
cg_idx=str2double(cg.RowLabels);


sort=idx(cg_idx);

input1=input(sort==sources(1),:);
input2=input(sort==sources(2),:);

cg1=clustergram(input1,'Cluster','column','Symmetric','false','Colormap','parula','ColumnLabels',targets);
cg2=clustergram(input2,'Cluster','column','Symmetric','false','Colormap','parula','ColumnLabels',targets);

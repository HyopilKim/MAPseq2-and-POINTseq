function drawCluster(Bnorm,targets,normhow)

%If you want to normalize cells by maximum projection, normhow=1, by total
%projection, normhow=0

if normhow==0
    clustergram(Bnorm,'Standardize','Row','ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')
else
    clustergram(Bnorm./repmat(max(Bnorm,[],2),1,size(Bnorm,2)),'Standardize',0,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')
end
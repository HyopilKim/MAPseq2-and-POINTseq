%% function to calculate the conditional projection probabilities
function [conditionalp,rowlabels,columnlables,matrix_binary,numerator]=conditionalprobabilities2(matrix,minneurons,cutoff,labels,plot)

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>minneurons);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
columnlables=labels(sum(matrix_binary_tmp,1)>minneurons);



rowlabels=columnlables;
j=1;
for i=1:length(columnlables)
    conditionalp(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1)./sum(matrix_binary(:,i)==1);
    numerator(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1);
    j=j+1;
end

if plot==1
    figure;imagesc(conditionalp');colorbar;
    ax = gca;
    ax.YTick=[1:length(columnlables)];
    ax.YTickLabel=columnlables;
    ax.XTick=[1:length(rowlabels)];
    ax.XTickLabel=rowlabels;
    ax.XTickLabelRotation=45;
    ylabel('P(B|A)');
    xlabel('Area A')
    
% let's divide P(B|A) with P(B) to check over or under representation
    divbyPb=conditionalp./sum(matrix_binary);
    figure;imagesc(divbyPb');colorbar;
    ax = gca;
    ax.YTick=[1:length(columnlables)];
    ax.YTickLabel=columnlables;
    ax.XTick=[1:length(rowlabels)];
    ax.XTickLabel=rowlabels;
    ax.XTickLabelRotation=45;
    ylabel('P(B|A)/P(B)')
    xlabel('Area A')
    
    
    
    % make a clustergram
    clustergram(conditionalp',...
    'Symmetric',0,'Standardize',3,'Cluster',3,'ColumnLabels',rowlabels,...
    'RowLabels',columnlables,'Colormap','parula',...
    'RowPDist','euclidean','ColumnPDist','euclidean');
end
findfigs;

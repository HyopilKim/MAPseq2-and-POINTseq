%% function to calculate the conditional projection probabilities
function [conditionalp,rowlabels,columnlables,matrix_binary,numerator,signif,coE,coO]=conditionalprobabilities3(matrix,minneurons,cutoff,labels,N_t,plot,label)

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>minneurons);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
columnlables=labels(sum(matrix_binary_tmp,1)>minneurons);

coP=zeros(length(columnlables));
coO=zeros(length(columnlables));
comb=nchoosek(1:length(columnlables),2);
P=sum(matrix_binary,1)/N_t;
for i=1:size(comb,1)
    r=comb(i,1);c=comb(i,2);
    coP(r,c)=P(r)*P(c);
    coO(r,c)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
end
coE=coP*N_t;



% Calculate the original p-value (two-tailed) for the observed outcome
pvalue=ones(length(columnlables));

for k=1:nchoosek(length(columnlables),2)
    r=comb(k,1);c=comb(k,2);
    if coO(r,c) > coE(r,c)
        pvalue(r,c) = 1 - binocdf(coO(r,c) - 1, N_t, coP(r,c));
    else
        pvalue(r,c) = binocdf(coO(r,c), N_t, coP(r,c));
    end
end


% Apply Bonferroni correction
p_adj = pvalue * nchoosek(length(columnlables),2);
signif=strings(size(p_adj));
signif(p_adj>=0.05)="";
signif(p_adj<0.05 & p_adj>=0.01)="*";
signif(p_adj<0.01 & p_adj>=0.001)="**";
signif(p_adj<0.001 & p_adj>=0.0001)="***";
signif(p_adj<0.0001)="****";
signif(p_adj==1)="";
for i=1:length(columnlables)
    signif(i,i)="";
end

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
    



signif=signif';
signif=signif(2:end,1:end-1);

triladj=tril(ones(length(columnlables)))';

% let's divide P(B|A) with P(B) to check over or under representation
    divbyPb=tril(N_t*conditionalp./sum(matrix_binary),-1);
    divbyPb=divbyPb+triladj;
    divbyPb=divbyPb(2:end,1:end-1);
    figure;imagesc(divbyPb);clim([0.0625 16]);colormap(redblue);colorbar;
    set(gca,'ColorScale','log');
    ax = gca;
    ax.YTick=1:(length(columnlables)-1);
    ax.YTickLabel=columnlables(2:end);
    ax.XTick=1:(length(rowlabels)-1);
    ax.XTickLabel=rowlabels(1:end-1);
    ax.XTickLabelRotation=45;
    [rows, cols] = size(divbyPb);
for i = 1:rows
   for j = 1:cols
       if j>i
       text(j, i, "", 'HorizontalAlignment', 'center');
       else
        if label==1
            text(j, i, num2str(divbyPb(i,j),'%.2f'), 'HorizontalAlignment', 'center');
            text(j, i-0.3, signif(i,j), 'HorizontalAlignment', 'center');
        else
            text(j, i, signif(i,j), 'HorizontalAlignment', 'center');
        end
       end
    end
end
    ylabel('P(A ∩ B)/P(A)P(B)');
    xlabel('Area A');
ax = gca;
ax.Box = 'off';       
ax.TickLength = [0 0]; 
    
    
    % make a clustergram
    clustergram(conditionalp',...
    'Symmetric',0,'Standardize',3,'Cluster',3,'ColumnLabels',rowlabels,...
    'RowLabels',columnlables,'Colormap','parula',...
    'RowPDist','euclidean','ColumnPDist','euclidean');
end

findfigs;


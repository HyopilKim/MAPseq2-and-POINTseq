
function conditionalP(matrix,targetlabels,cutoff,min_region,min_cts,dolabel,doN_t)


n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
Plabels=targetlabels(sum(matrix_binary_tmp,1)>=min_region);


%find N_t counts
if exist('doN_t','var')
    [N_obs, k] = size(matrix_binary);
    N_a = sum(matrix_binary, 1);  % count 1s in each dimension

% Initialize polynomial coefficients array
    coeffs = zeros(1, k+1);
   
% First coefficient (highest degree)
    coeffs(1) = N_obs - sum(N_a);
   
% Calculate other coefficients using nchoosek
    for i = 2:k
   % (-1)^(i-1) gives alternating signs
        coeffs(i) = -(-1)^(i-1) * sum(prod(nchoosek(N_a,i),2));
    end
   
% Last coefficient (constant term)
    coeffs(end) = -(-1)^k * prod(N_a);
   
% Find roots
    r = roots(coeffs);
   
% Take the smallest positive real root as N_total
    N_t = round(max(r));
else
    N_t=size(matrix_binary,2);
end


coP=zeros(length(Plabels));
coO=zeros(length(Plabels));
comb=nchoosek(1:length(Plabels),2);
P=sum(matrix_binary,1)/N_t;
for i=1:size(comb,1)
    r=comb(i,1);c=comb(i,2);
    coP(r,c)=P(r)*P(c);
    coO(r,c)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
    sum(sum(matrix_binary(:,[r c]),2)>0);
    
    sum(matrix_binary,1)
end
coE=coP*N_t;






% Calculate the original p-value (two-tailed) for the observed outcome
pvalueP=ones(length(Plabels));

for k=1:nchoosek(length(Plabels),2)
    r=comb(k,1);c=comb(k,2);
    if coO(r,c) > coE(r,c)
        pvalueP(r,c) = 1 - binocdf(coO(r,c) - 1, N_t, coP(r,c));
    else
        pvalueP(r,c) = binocdf(coO(r,c), N_t, coP(r,c));
    end
end


% Apply Bonferroni correction
p_adjP = pvalueP * nchoosek(length(Plabels),2);
pre_signifP=strings(size(p_adjP));
pre_signifP(p_adjP>=0.05)="NS";
pre_signifP(p_adjP<0.05 & p_adjP>=0.01)="*";
pre_signifP(p_adjP<0.01 & p_adjP>=0.001)="**";
pre_signifP(p_adjP<0.001 & p_adjP>=0.0001)="***";
pre_signifP(p_adjP<0.0001)="****";
pre_signifP(p_adjP==1)="";


for i=1:length(Plabels)
    for j=1:length(Plabels)
        pre_signifP(max(coE(j,i), coO(j,i))<min_cts)="";
    end
    pre_signifP(i,i)="";
end


[p_adjP_idxr, p_adjP_idxc]=find((p_adjP<0.05)==1);
p_adjP_idx=[p_adjP_idxr p_adjP_idxc];


j=1;


conditionalp=zeros(length(Plabels));
numerator=zeros(length(Plabels));
for i=1:length(Plabels)
    conditionalp(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1)./sum(matrix_binary(:,i)==1);
    numerator(j,:)=sum(matrix_binary(matrix_binary(:,i)==1,:),1);
    j=j+1;
end


figure;imagesc(conditionalp');clim([0 1]);colorbar;
ax = gca;
ax.YTick=1:length(Plabels);
ax.YTickLabel=Plabels;
ax.XTick=1:length(Plabels);
ax.XTickLabel=Plabels;
ax.XTickLabelRotation=45;
ylabel('P(B|A)');
xlabel('Area A')
    



signifP=pre_signifP';
signifP=signifP(2:end,1:end-1);

triladj=tril(ones(length(Plabels)))';

% let's divide P(B|A) with P(B) to check over or under representation
divbyPb=tril(N_t*conditionalp./sum(matrix_binary),-1);
divbyPb=divbyPb+triladj;
divbyPb=divbyPb(2:end,1:end-1);
    
figure;imagesc(divbyPb);clim([0.0625 16]);colormap(redblue);colorbar;
set(gca,'ColorScale','log');
ax = gca;
ax.YTick=1:(length(Plabels)-1);
ax.YTickLabel=Plabels(2:end);
ax.XTick=1:(length(Plabels)-1);
ax.XTickLabel=Plabels(1:end-1);
ax.XTickLabelRotation=45;
[rows, cols] = size(divbyPb);
for i = 1:rows
   for j = 1:cols
       if j>i
       text(j, i, "", 'HorizontalAlignment', 'center');
       else
        if dolabel==1
            text(j, i, num2str(divbyPb(i,j),'%.2f'), 'HorizontalAlignment', 'center');
            text(j, i-0.3, signifP(i,j), 'HorizontalAlignment', 'center');
        else
            text(j, i, signifP(i,j), 'HorizontalAlignment', 'center');
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
conditionalp_t=conditionalp';
    
cg = clustergram(conditionalp_t);
rowOrder = cg.RowLabels; 
rowIdx=zeros(1,length(rowOrder));
for i=1:length(rowOrder)
    rowIdx(i)=str2num(rowOrder{i});
end
colOrder = cg.ColumnLabels;
colIdx=zeros(1,length(colOrder));
for i=1:length(colOrder)
    colIdx(i)=str2num(colOrder{i});
end
    
    
reorderedP = conditionalp_t(rowIdx, colIdx);
figure;imagesc(reorderedP);clim([0 1]);colorbar;
ax = gca;
ax.YTick=1:length(Plabels);
ax.YTickLabel=Plabels(rowIdx);
ax.XTick=1:length(Plabels);
ax.XTickLabel=Plabels(colIdx);
ax.XTickLabelRotation=45;
ylabel('P(B|A)');
xlabel('Area A')

findfigs;

save('contionalP.mat','conditionalp','matrix_binary','N_t','signifP','coE','coO','p_adjP');
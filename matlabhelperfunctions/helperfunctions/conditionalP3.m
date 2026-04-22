
function conditionalP3(matrix,targetlabels,cutoff,min_region,min_cts,pscale,uscale,nscale,dscale,cscale,dolabel,doN_t)

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);

matrix1=matrix.*matrix_binary_tmp;
matrix2=matrix1(:,sum(matrix_binary_tmp,1)>=min_region);
matrix3=matrix2(sum(matrix_binary_tmp2,2)~=0,:);

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
    N_t=size(matrix_binary,1);
end


coP=zeros(length(Plabels));
coO=zeros(length(Plabels));
comb=nchoosek(1:length(Plabels),2);
P=sum(matrix_binary,1)/N_t;
for i=1:size(comb,1)
    r=comb(i,1);c=comb(i,2);
    coP(r,c)=P(r)*P(c);
    coP(c,r)=P(r)*P(c);
    coO(r,c)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
    coO(c,r)=sum((matrix_binary(:,r)>0) & (matrix_binary(:,c)>0));
end
coE=coP*N_t;



% Calculate the original p-value (two-tailed) for the observed outcome
pvalueP=ones(length(Plabels));

for k=1:nchoosek(length(Plabels),2)
    r=comb(k,1);c=comb(k,2);
    if coO(r,c) > coE(r,c)
        pvalueP(r,c) = 1 - binocdf(coO(r,c) - 1, N_t, coP(r,c));
        pvalueP(c,r) = 1 - binocdf(coO(r,c) - 1, N_t, coP(r,c));
    else
        pvalueP(r,c) = binocdf(coO(r,c), N_t, coP(r,c));
        pvalueP(c,r) = binocdf(coO(r,c), N_t, coP(r,c));
    end
end



% Apply Bonferroni correction
p_adjP = pvalueP * nchoosek(length(Plabels),2);
signifP=strings(size(p_adjP));
signifP(p_adjP>=0.05)="NS";
signifP(p_adjP<0.05 & p_adjP>=0.01)="*";
signifP(p_adjP<0.01 & p_adjP>=0.001)="**";
signifP(p_adjP<0.001 & p_adjP>=0.0001)="***";
signifP(p_adjP<0.0001)="****";
signifP(p_adjP==1)="";


for i=1:length(Plabels)
    for j=1:length(Plabels)
        if max(coE(j,i), coO(j,i))<min_cts
            signifP(j,i)="Few";
        end
    end
    signifP(i,i)="";
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

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(conditionalp');clim([0 pscale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(conditionalp')
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;


conditionalp_u=zeros(length(Plabels));
for i=1:length(Plabels)
    for j=1:length(Plabels)
        conditionalp_u(i,j)=length(matrix_binary(matrix_binary(:,i)==1 & matrix_binary(:,j)==1))/length(matrix_binary(matrix_binary(:,i)==1 | matrix_binary(:,j)==1));
    end
end
    

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(conditionalp_u);clim([0 uscale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(conditionalp_u)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
imagesc(numerator/N_t);clim([0 nscale]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(numerator)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

% let's divide P(B|A) with P(B) to check over or under representation
divbyPb=N_t*conditionalp./sum(matrix_binary);

figure;imagesc(divbyPb);clim(dscale);colormap(redblue);
set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
set(gca,'ColorScale','log');
ax = gca;
ax.YTick=[];
ax.XTick=[];
[rows, cols] = size(divbyPb);
for i = 1:rows
   for j = 1:cols
        if dolabel==1
        text(j, i, num2str(divbyPb(i,j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',10);
        text(j, i-0.3, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        else
        text(j, i, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        end
    end
end
hold on;
for i = 1:length(divbyPb)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

%correlation of input matrix    
[pre_corrMatrix, pMatrix] = corr(matrix3, 'Type', 'Spearman');
corrMatrix=abs(pre_corrMatrix);
corrSig=strings(size(pMatrix));
corrSig(pMatrix>=0.05)="NS";
corrSig(pMatrix<0.05 & pMatrix>=0.01)="*";
corrSig(pMatrix<0.01 & pMatrix>=0.001)="**";
corrSig(pMatrix<0.001 & pMatrix>=0.0001)="***";
corrSig(pMatrix<0.0001)="****";

figure;imagesc(corrMatrix);clim(cscale);colormap(redblue);
set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
[rows, cols] = size(corrMatrix);
for i = 1:rows
   for j = 1:cols
        if dolabel==1
        text(j, i, num2str(corrMatrix(i,j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',10);
        text(j, i-0.3, corrSig(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        else
        text(j, i, corrSig(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        end
    end
end
hold on;
for i = 1:length(corrMatrix)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;



%correlation of input matrix (nonzero)    
[corrMatrix, pMatrix] = nonzero_pairwise_corr(matrix3, 'Pearson');
corrSig=strings(size(pMatrix));
corrSig(pMatrix>=0.05)="NS";
corrSig(pMatrix<0.05 & pMatrix>=0.01)="*";
corrSig(pMatrix<0.01 & pMatrix>=0.001)="**";
corrSig(pMatrix<0.001 & pMatrix>=0.0001)="***";
corrSig(pMatrix<0.0001)="****";

figure;imagesc(corrMatrix);clim(cscale);colormap(redblue);
set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);
ax = gca;
ax.YTick=[];
ax.XTick=[];
[rows, cols] = size(corrMatrix);
for i = 1:rows
   for j = 1:cols
        if dolabel==1
        text(j, i, num2str(corrMatrix(i,j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',10);
        text(j, i-0.3, corrSig(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        else
        text(j, i, corrSig(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        end
    end
end
hold on;
for i = 1:length(corrMatrix)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;


findfigs;

save('contionalP.mat','conditionalp','conditionalp_u','matrix_binary','N_t','signifP','coE','coO','p_adjP','Plabels','numerator');
%% Return observed counts or estimated 
function [obs,esti,p_adj,signif,relative_obs]=motif(matrix_binary,N_t,maxtar,targets,minneuron)

%counts for observed motifs
n=size(matrix_binary,2);
obs=[];
k=1;
for i=2:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)    
    obs(k)=sum(min(matrix_binary(:,comb(j,:)),[],2)==1 & sum(matrix_binary,2)==i);
    k=k+1;
    end
end

%counts for estimated motifs
esti=[];
k=1;
P=sum(matrix_binary)/N_t;
noP=ones(1,size(matrix_binary,2))-P;
for i=2:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)
        esti(k)=N_t*prod(P(comb(j,:)))*prod(noP(setdiff(1:n, comb(j))));
        k=k+1;
    end
end

precolnum=length(obs);

filt=[];
if exist('minneuron','var')
    for i=1:length(obs)
        if (max(obs(i),esti(i))>minneuron)
            filt=[filt i];
        end
    end
    obs=obs(filt);
    esti=esti(filt);
end
        
% Calculate the original p-value (two-tailed) for the observed outcome
pvalue=zeros(1,length(obs));
for k=1:length(obs)
    if obs(k) > esti(k)
        pvalue(k) = 1 - binocdf(obs(k) - 1, N_t, esti(k)/N_t);
    else
        pvalue(k) = binocdf(obs(k), N_t, esti(k)/N_t);
    end
end


% Apply Bonferroni correction
p_adj = pvalue * length(obs);
signif=strings(size(p_adj));
signif(p_adj>=0.05)="NS";
signif(p_adj<0.05 & p_adj>=0.01)="*";
signif(p_adj<0.01 & p_adj>=0.001)="**";
signif(p_adj<0.001 & p_adj>=0.0001)="***";
signif(p_adj<0.0001)="****";

%plot
plot_target=zeros(n,precolnum);
k=1;
for i=2:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)
        plot_target(comb(j,:),k)=1;
        k=k+1;
    end
end

if exist('minneuron','var')
    plot_target=plot_target(1:n,filt);
end

relative_obs=obs./esti;

str_esti = cell(1, numel(esti));
for i = 1:numel(esti)
    if esti(i) < 0.05
        str_esti{i} = sprintf('%.0e', esti(i));  
    else
        str_esti{i} = sprintf('%.1f', esti(i));  
    end
end
str_esti = string(str_esti);  

sig_counts=string(signif)+newline+string(obs)+"/"+newline+str_esti;

sig_counts_big=sig_counts;
sig_counts_big((relative_obs<0.05)|(relative_obs>50))="";

sig_counts_small=sig_counts;
sig_counts_small(relative_obs>=0.05)="";

sig_counts_super=sig_counts;
sig_counts_super(relative_obs<=50)="";

overunder_colors = zeros(length(relative_obs), 3); 
overunder_colors((relative_obs > 1) & (p_adj<0.05),:) = repmat([1 0 0],size(overunder_colors((relative_obs > 1) & (p_adj<0.05),:),1),1);
overunder_colors((relative_obs > 1) & (p_adj>=0.05),:) = repmat([0.5 0.5 0.5],size(overunder_colors((relative_obs > 1) & (p_adj>=0.05),:),1),1);
overunder_colors((relative_obs < 1) & (p_adj<0.05),:) = repmat([0 0 1],size(overunder_colors((relative_obs < 1) & (p_adj<0.05),:),1),1);
overunder_colors((relative_obs < 1) & (p_adj>=0.05),:) = repmat([0.5 0.5 0.5],size(overunder_colors((relative_obs < 1) & (p_adj>=0.05),:),1),1);






figure;set(gcf,'Units', 'normalized', 'Position', [0, 0.05, (0.9), 0.6]);
ax1=subplot(24,1,1:17);
bar(relative_obs,'FaceColor','flat','CData',overunder_colors);
set(gca,'XTick',[],'FontSize',10,'YScale','log','YLim',[0.05 50],'YTick', [0.1 1 10], 'YTickLabel', {'0.1', '1', '10'});
text(1:length(relative_obs), relative_obs*2, sig_counts_big,'FontWeight','bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'top','FontSize',7);
text(1:length(relative_obs), repmat(0.1,1,length(relative_obs)),sig_counts_small,'FontWeight','bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'top','FontSize',7);
text(1:length(relative_obs), repmat(40,1,length(relative_obs)),sig_counts_super,'FontWeight','bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'top','FontSize',7);
ax2=subplot(24,1,18:24);
imagesc(plot_target);set(gcf,'Colormap',colormap(flipud(gray)));
set(gca,'XTick',[]);
if exist('targets','var')
    set(gca,'YTick',1:length(targets), 'YTickLabels', targets);
else
    set(gca,'YTick',[]);
end
linkaxes([ax1 ax2],'x');
findfigs;


n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

figure;set(gcf,'Units', 'normalized', 'Position', [0, 0.05, (0.9*length(obs)/56), 0.6]);
ax1=subplot(24,1,1:17);
bar(1:length(obs),[obs;esti]);
set(gca,'XTick',[],'FontSize',10);
ax2=subplot(24,1,18:24);
imagesc(plot_target);set(gcf,'Colormap',colormap(flipud(gray)));
set(gca,'XTick',[]);
if exist('targets','var')
    set(gca,'YTick',1:length(targets), 'YTickLabels', targets);
else
    set(gca,'YTick',[]);
end
linkaxes([ax1 ax2],'x');
findfigs;



plot_target_color=ones(size(plot_target));
for i=1:size(plot_target_color,2)
    plot_target_color(plot_target(:,i)==1,i)=relative_obs(i);
end

signif_filt=signif(p_adj<0.05);

relative_filt=relative_obs(p_adj<0.05);

figure;set(gcf,'Units', 'normalized', 'Position', [0, 0.05, (sum(p_adj<0.05)/27)*(0.9*length(obs)/56), 0.6]);
ax1=subplot(30,1,1:17);
obs_filt=obs(p_adj<0.05);
esti_filt=esti(p_adj<0.05);
bar(1:length(obs_filt),[obs_filt;esti_filt]);
set(gca,'XTick',[],'YLim',[0 55],'FontSize',12);
ax2=subplot(30,1,18:30);
plot_target_filt=plot_target_color(:,p_adj<0.05);
imagesc(plot_target_filt);clim([1/64 64]);colormap(redblue);set(gca,'ColorScale','log');
for i=1:length(targets)
    for j=1:sum(p_adj<0.05)
        if plot_target_filt(i,j)==1
            text(j, i, "", 'HorizontalAlignment', 'center');
            text(j, i-0.3, "", 'HorizontalAlignment', 'center');
        else
            text(j, i, num2str(plot_target_filt(i,j),'%.2f'), 'HorizontalAlignment', 'center');
            text(j, i-0.3, signif_filt(j), 'HorizontalAlignment', 'center');
        end
    end
end
if exist('targets','var')
    set(gca,'YTick',1:length(targets), 'YTickLabels', targets);
else
    set(gca,'YTick',[]);
end
linkaxes([ax1 ax2],'x');
set(gca,'XTick',[]);
findfigs;

figure;set(gcf,'Units', 'normalized', 'Position', [0, 0.05, (sum(p_adj<0.05)/27)*(0.9*length(obs)/56), 0.6]);
ax1=subplot(30,1,1:17);
obs_filt=obs(p_adj<0.05);
esti_filt=esti(p_adj<0.05);
bar(1:length(obs_filt),[obs_filt;esti_filt]);
set(gca,'XTick',[],'YLim',[0 350],'FontSize',12);
ax2=subplot(30,1,18:30);
plot_target_filt=plot_target_color(:,p_adj<0.05);
imagesc(plot_target_filt);clim([1/64 64]);colormap(redblue);set(gca,'ColorScale','log');
for j=1:sum(p_adj<0.05)
    text(j, length(targets)+1, num2str(relative_filt(j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',14);
    text(j, length(targets)+1-0.3, signif_filt(j), 'HorizontalAlignment', 'center','FontSize',14);
end
if exist('targets','var')
    set(gca,'YTick',1:length(targets), 'YTickLabels', targets,'YAxisLocation', 'right');
else
    set(gca,'YTick',[]);
end
linkaxes([ax1 ax2],'x');
set(gca,'XTick',[]);
findfigs;
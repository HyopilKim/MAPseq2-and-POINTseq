function motif2(matrix,targetlabels,cutoff,cut_relativeobs,min_region,min_cts,p_adjP,maxtar,N_t,doeachneuron)

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

[p_adjP_idxr, p_adjP_idxc]=find((p_adjP<0.05)==1);
p_adjP_idx=[p_adjP_idxr p_adjP_idxc];

%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;

matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
Mlabels=targetlabels(sum(matrix_binary_tmp,1)>=min_region);

matrix_filt1=matrix.*(matrix>cutoff);
matrix_filt=matrix_filt1(sum(matrix_binary_tmp2,2)~=0,:);


n=size(matrix_binary,2);
pre_neurons_mtf={};
pre_norm_neurons_mtf={};
prepre_obs=[];
k=1;

for i=2:maxtar
    comb=nchoosek(1:length(targetlabels),i);
    for j=1:size(comb,1)    
    except=ones(1,size(matrix_binary,2));
    except(comb(j,:))=0;
    prepre_neurons_mtf_tmp=matrix_filt(sum(matrix_binary.*except,2)==0,:);
    prepre_neurons_mtf{k}=prepre_neurons_mtf_tmp;
    prepre_neurons_mtfonly{k}=matrix_filt(min(matrix_binary(:,comb(j,:)),[],2)==1 & sum(matrix_binary,2)==i,:);
    prepre_norm_neurons_mtf{k}=prepre_neurons_mtf_tmp./max(prepre_neurons_mtf_tmp,[],2);
    prepre_norm_neurons_mtfonly{k}=prepre_neurons_mtfonly{k}./max(prepre_neurons_mtfonly{k},[],2);
    prepre_obs(k)=sum(min(matrix_binary(:,comb(j,:)),[],2)==1 & sum(matrix_binary,2)==i);
    k=k+1;
    end
end


prepre_esti=[];
k=1;
P=sum(matrix_binary)/N_t;
noP=ones(1,size(matrix_binary,2))-P;
for i=2:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)
        prepre_esti(k)=N_t*prod(P(comb(j,:)))*prod(noP(setdiff(1:n, comb(j))));
        k=k+1;
    end
end

filt=[];
if exist('min_cts','var')
    for i=1:length(prepre_obs)
        if (max(prepre_obs(i),prepre_esti(i))>=min_cts)
            filt=[filt i];
        end
    end
    pre_obs=prepre_obs(filt);
    pre_esti=prepre_esti(filt);
    pre_relative_obs=pre_obs./pre_esti;
    pre_neurons_mtf=prepre_neurons_mtf(filt);
    pre_norm_neurons_mtf=prepre_norm_neurons_mtf(filt);
    pre_neurons_mtfonly=prepre_neurons_mtfonly(filt);
    pre_norm_neurons_mtfonly=prepre_norm_neurons_mtfonly(filt);
end

relative_obs=pre_relative_obs(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
obs=pre_obs(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
esti=pre_esti(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));

neurons_mtf=pre_neurons_mtf(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
norm_neurons_mtf=pre_norm_neurons_mtf(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
neurons_mtfonly=pre_neurons_mtfonly(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
norm_neurons_mtfonly=pre_norm_neurons_mtfonly(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));


% Calculate the original p-value (two-tailed) for the observed outcome
pre_p = pre_esti / N_t;

% 결과 벡터 초기화
pre_pvalueM = nan(size(pre_obs));

% 각 원소에 대해 binomial two-tailed p-value 계산
for i = 1:length(pre_obs)
    k_obs = pre_obs(i);
    p_null = pre_p(i);

    x = 0:N_t;
    prob = binopdf(x, N_t, p_null);
    p_k = binopdf(k_obs, N_t, p_null);

    % two-tailed p-value
    pre_pvalueM(i) = sum(prob(prob <= p_k));
end

pvalueM=pre_pvalueM(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));

% Apply Bonferroni correction
pre_p_adjM = pre_pvalueM * length(obs);
p_adjM = pvalueM * length(obs);
signifM=strings(size(p_adjM));
signifM(p_adjM>=0.05)="NS";
signifM(p_adjM<0.05 & p_adjM>=0.01)="*";
signifM(p_adjM<0.01 & p_adjM>=0.001)="**";
signifM(p_adjM<0.001 & p_adjM>=0.0001)="***";
signifM(p_adjM<0.0001)="****";


%plot
prepre_plot_target=zeros(length(targetlabels),length(prepre_obs));
k=1;
for i=2:maxtar
    comb=nchoosek(1:length(targetlabels),i);
    for j=1:size(comb,1)
        prepre_plot_target(comb(j,:),k)=1;
        k=k+1;
    end
end

pre_plot_target=prepre_plot_target(:,filt);
pre_plot_target_color=repmat(pre_relative_obs,n,1);
pre_plot_target_color(pre_plot_target==0)=1;

plot_target=pre_plot_target(:,pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
plot_target_color=pre_plot_target_color(:,pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));


mark_targets_all=repmat(Mlabels',1,size(plot_target,2));
mark_targets=strings(size(mark_targets_all));
mark_targets(logical(plot_target))=mark_targets_all(logical(plot_target));
mark_targets_filt=mark_targets(:,p_adjM<0.05);

Pconnected_motif_asidx={};
Pconnected_motif_aslabel={};
coRegions="";
for i=1:size(p_adjP_idx,1)
    selectbyP=(p_adjM<0.05)&(plot_target(p_adjP_idx(i,1),:)==1 & plot_target(p_adjP_idx(i,2),:)==1);
    asidx=plot_target(:,selectbyP);
    Pconnected_motif_asidx{i}=asidx;
    coRegions(i)=strjoin(Mlabels(p_adjP_idx(i,:)),'-');
    aslabel="";
    for j=1:size(asidx,2)
    aslabel(j)=strjoin(Mlabels(logical(asidx(:,j))),'-');
    end
    Pconnected_motif_aslabel{i}=aslabel;
end
   

str_esti = cell(1, numel(esti));
for i = 1:numel(esti)
    if esti(i) < 0.05
        str_esti{i} = sprintf('%.0e', esti(i));  
    else
        str_esti{i} = sprintf('%.1f', esti(i));  
    end
end
str_esti = string(str_esti);  

sig_counts=string(signifM)+newline+string(obs)+"/"+newline+str_esti;

sig_counts_big=sig_counts;
sig_counts_big((relative_obs<0.05)|(relative_obs>50))="";

sig_counts_small=sig_counts;
sig_counts_small(relative_obs>=0.05)="";

sig_counts_super=sig_counts;
sig_counts_super(relative_obs<=50)="";

relative_filt=round(relative_obs(:,p_adjM<0.05),2);
relative_filt_str=string(relative_filt);
relative_filt_str(relative_filt>100)=">100";

signifM_filt=signifM(p_adjM<0.05);

neurons_mtf_filt=neurons_mtf(p_adjM<0.05);
norm_neurons_mtf_filt=norm_neurons_mtf(p_adjM<0.05);
neurons_mtfonly_filt=neurons_mtfonly(p_adjM<0.05);
norm_neurons_mtfonly_filt=norm_neurons_mtfonly(p_adjM<0.05);


%motif plot 
figure;set(gcf,'Units', 'normalized', 'Position', [0, 0,sum(p_adjM<0.05)*0.04, 0.5]);
ax1=subplot(30,1,1:16);
obs_filt=obs(p_adjM<0.05);
esti_filt=esti(p_adjM<0.05);
bar(1:length(obs_filt),[obs_filt;esti_filt]);
motif_ylim=max([obs_filt esti_filt])*1.1;
set(gca,'XTick',[],'YLim',[0 motif_ylim],'FontSize',12);
for j=1:sum(p_adjM<0.05)
    text(j, max(obs_filt(j),esti_filt(j))+motif_ylim/30, relative_filt_str(j), 'HorizontalAlignment', 'center','FontSize',13);
    text(j, max(obs_filt(j),esti_filt(j))+motif_ylim/15, signifM_filt(j), 'HorizontalAlignment', 'center','FontSize',13,'FontWeight','bold');
end
ax2=subplot(30,1,17:25);
plot_target_filt=plot_target_color(:,p_adjM<0.05);
imagesc(plot_target_filt);clim([1/64 64]);colormap(redblue);set(gca,'ColorScale','log');
for i=1:length(Mlabels)
    for j=1:sum(p_adjM<0.05)
        text(j, 1+length(Mlabels), string(j), 'HorizontalAlignment', 'center','FontSize',15);
        text(j, i, mark_targets_filt(i,j), 'HorizontalAlignment', 'center','FontSize',11,'FontWeight','bold');
    end
end
linkaxes([ax1 ax2],'x');
set(gca,'XTick',[],'YTick',[]);
findfigs;

%volcano plot
if min(pre_relative_obs)<1
    relobsnonzero=pre_relative_obs;
    relobsnonzero(pre_relative_obs==0)=max(pre_relative_obs);
    minrelobs=min(relobsnonzero);
else
    minrelobs=0;
end

relative_obs_volcano=log2(pre_relative_obs+minrelobs);
p_volcano=-log10(pre_p_adjM+min(pre_p_adjM(pre_p_adjM~=0)));

pre_text_x=relative_obs_volcano(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
pre_text_y=p_volcano(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));

figure;scatter(relative_obs_volcano,p_volcano,50,'b', 'filled');hold on;
yline(-log10(0.05+min(pre_p_adjM)),'--r','LineWidth',1);
xline(log2(1),'--r','LineWidth',1);
text(pre_text_x(p_adjM<0.05),pre_text_y(p_adjM<0.05),string(1:sum(p_adjM<0.05)),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',15);
hold off;
if sum(isinf([log2(relative_obs) log10(pre_p_adjM)]))>0 
    set(gca,'FontSize',13,'XLim',[min(relative_obs_volcano) max(relative_obs_volcano)],'YLim',[min(p_volcano) max(p_volcano)]);
else
    set(gca,'FontSize',13);
end




%%individual neurons of each motif plot (normalized)
mark_mtf="";
for i=1:length(mark_targets_filt)
    mark_tmp=mark_targets_filt(:,i);
    mark_tmp=mark_tmp(mark_tmp~="");
    mark_mtf(i)=strjoin(mark_tmp,"-");
end

if exist('doeachneuron','var')
    %as line plot
    %figure;
    %for i=1:length(norm_neurons_mtf_filt)
       % plot(1:length(targetlabels),norm_neurons_mtf_filt{i},'-o','LineWidth',2);set(gca,'XTickLabel',targetlabels);
        %findfigs;exportgraphics(gcf, sprintf('eachneuron_mtf_line_%d.png', i));
   % end

   % figure;
   % for i=1:length(norm_neurons_mtfonly_filt)
      %  if ~isempty(norm_neurons_mtfonly_filt{i})
           % plot(1:length(targetlabels),norm_neurons_mtfonly_filt{i},'-o','LineWidth',2);set(gca,'XTickLabel',targetlabels);
           % findfigs;exportgraphics(gcf, sprintf('eachneuron_mtfonly_line_%d.png', i));
       % end
   % end

    %as matrix
   % targetfilt=plot_target(:,p_adjM<0.05);

   % figure;
    %for i=1:length(norm_neurons_mtf_filt)
       % yfactor=(length(norm_neurons_mtf_filt{i})/size(matrix_filt,1))*0.8;
       % set(gcf,'Units', 'normalized', 'Position', [0.1, 0.05, 0.45, 0.05+yfactor]);
       % imagesc(norm_neurons_mtf_filt{i});set(gca,'XTickLabel',targetlabels);colormap('parula');
      %  findfigs;exportgraphics(gcf, sprintf('eachneuron_mtf_matrix_%d.png', i));
    %end

 %   figure;
 %   for i=1:length(norm_neurons_mtfonly_filt)
 %       if ~isempty(norm_neurons_mtfonly_filt{i})
 %           yfactor=(length(norm_neurons_mtfonly_filt{i})/size(matrix_filt,1))*10;
 %           set(gcf,'Units', 'normalized', 'Position', [0.1, 0.05, 0.45, 0.05+yfactor]);
  %          imagesc(norm_neurons_mtfonly_filt{i});set(gca,'XTickLabel',targetlabels);colormap('parula');
   %         title(mark_mtf(i),'FontSize',18)
    %        ax = gca;ax.YAxis.FontSize = 14;ax.XAxis.FontSize = 16;ax.PositionConstraint = 'outerposition';
     %       row_count = size(norm_neurons_mtfonly_filt{i}, 1);
      %      if row_count <= 5
       %         yticks(0:2:row_count);  
        %    elseif row_count < 30
         %       yticks(0:5:row_count);
          %  elseif row_count < 200
           %     yticks(0:20:row_count); 
            %elseif row_count < 500
            %    yticks(0:50:row_count);  
            %elseif row_count < 5000
            %    yticks(0:500:row_count);
            %else
            %    yticks(linspace(0, row_count, 10)); 
            %end
            %findfigs;pause(1);exportgraphics(gcf, sprintf('eachneuron_mtfonly_matrix_%d.png', i));
        %end
    %end

    figure;
    for i=1:length(norm_neurons_mtfonly_filt)
        if ~isempty(norm_neurons_mtfonly_filt{i})
            set(gcf,'Units', 'normalized', 'OuterPosition', [0, 0, 0.5, 0.3]);
            imagesc(norm_neurons_mtfonly_filt{i});set(gca,'XTickLabel',targetlabels);colormap('parula');
            title(string(i)+". "+mark_mtf(i),'FontSize',24)
            ax = gca;ax.YAxis.FontSize = 18;ax.XAxis.FontSize = 20;ax.PositionConstraint = 'outerposition';
            findfigs;pause(1);exportgraphics(gcf, sprintf('eachneuron_mtfonly_matrix_%d.png', i));
        end
    end
end

findfigs;

writematrix(coRegions','coprojection and motif.xlsx','Sheet','Sheet1');
for i=1:length(Pconnected_motif_aslabel)
    writematrix(Pconnected_motif_aslabel{i},'coprojection and motif.xlsx','Sheet','Sheet1','Range',("B"+num2str(i)));
end

save('motif.mat','matrix_binary','obs','esti','p_adjM','signifM','relative_obs','Pconnected_motif_asidx','Pconnected_motif_aslabel','norm_neurons_mtf_filt','norm_neurons_mtfonly_filt');
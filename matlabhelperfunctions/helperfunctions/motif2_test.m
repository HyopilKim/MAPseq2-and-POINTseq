function motif2(matrix,targetlabels,cutoff,min_region,min_cts,p_adjP,maxtar,N_t_input,volcanoplotsize,volcanodotlabel,doeachneuron)


%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;
matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
Mlabels=targetlabels(sum(matrix_binary_tmp,1)>=min_region);


if N_t_input==0
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
   
% Take the largest positive real root as N_total
    N_t = round(max(r));
else
    N_t=size(matrix_binary,1);
end

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

[p_adjP_idxr, p_adjP_idxc]=find((p_adjP<0.05)==1);
p_adjP_idx=[p_adjP_idxr p_adjP_idxc];

%observed counts across all the motif combinations
n=size(matrix_binary,2);
pre_obs=[];
motif_idx={};
k=1;
for i=1:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)    
        motif_idx{k}=sum(matrix_binary(:,comb(j,:)),2)==length(comb(j,:));
        pre_obs(k)=sum(motif_idx{k});
        k=k+1;
    end
end

%expected counts across all the motif combinations
pre_exp=[];
k=1;
P=sum(matrix_binary)/N_t;
noP=ones(1,size(matrix_binary,2))-P;
for i=1:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)
        pre_exp(k)=N_t*prod(P(comb(j,:)))*prod(noP(setdiff(1:n, comb(j,:))));
        k=k+1;
    end
end

%Filter motifs >= min_cts for data reliability
filt=[];
if exist('min_cts','var')
    for i=1:length(pre_obs)
        if (max(pre_obs(i),pre_exp(i))>=min_cts)
            filt=[filt i];
        end
    end
    obs=pre_obs(filt);
    exp=pre_exp(filt);
    relative_obs=obs./exp;
end


% pvalue - binomial, one-tailed
pvalueM = nan(size(obs));
for i = 1:length(obs)
    if obs(i) > exp(i)
        pvalueM(i) = 1 - binocdf(obs(i) - 1, N_t, exp(i)/N_t);
    else
        pvalueM(i) = binocdf(obs(i), N_t, exp(i)/N_t);
    end
end

p_adjM = pvalueM * length(obs);% Apply Bonferroni correction
signifM=strings(size(p_adjM));
signifM(p_adjM>=0.05)="NS";
signifM(p_adjM<0.05 & p_adjM>=0.01)="*";
signifM(p_adjM<0.01 & p_adjM>=0.001)="**";
signifM(p_adjM<0.001 & p_adjM>=0.0001)="***";
signifM(p_adjM<0.0001)="****";

%plot targets
pre_plot_target=zeros(n,length(pre_obs));
k=1;
for i=1:maxtar
    comb=nchoosek(1:n,i);
    for j=1:size(comb,1)
        pre_plot_target(comb(j,:),k)=1;
        k=k+1;
    end
end
plot_target=pre_plot_target(:,filt);
plot_target_color=repmat(relative_obs,n,1);
plot_target_color(plot_target==0)=1;

mark_targets_all=repmat(Mlabels',1,size(plot_target,2));
mark_targets=strings(size(mark_targets_all));
mark_targets(logical(plot_target))=mark_targets_all(logical(plot_target));
mark_targets_filt=mark_targets(:,p_adjM<0.05);

%Motifs belongs to co-innervations
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
   
str_exp = cell(1, numel(exp));
for i = 1:numel(exp)
    if exp(i) < 0.05
        str_exp{i} = sprintf('%.0e', exp(i));  
    else
        str_exp{i} = sprintf('%.1f', exp(i));  
    end
end


relative_filt=round(relative_obs(:,p_adjM<0.05),2);
relative_filt_str=string(relative_filt);
relative_filt_str(relative_filt>100)=">100";




signifM_filt=signifM(p_adjM<0.05);

neurons_mtf_filt=neurons_mtf(p_adjM<0.05);
norm_neurons_mtf_filt=norm_neurons_mtf(p_adjM<0.05);
neurons_mtfonly_filt=neurons_mtfonly(p_adjM<0.05);
norm_neurons_mtfonly_filt=norm_neurons_mtfonly(p_adjM<0.05);

if sum(p_adjM<0.05)==0
    disp('no significant motif')
else
%motif plot 
figure;set(gcf,'Units', 'normalized', 'Position', [0, 0,sum(p_adjM<0.05)*0.03, 0.5]);
ax1=subplot(30,1,1:16);
obs_filt=obs(p_adjM<0.05);
exp_filt=exp(p_adjM<0.05);
bar(1:length(obs_filt),[obs_filt;exp_filt]);
motif_ylim=max([obs_filt exp_filt])*1.1;
set(gca,'XTick',[],'YLim',[0 motif_ylim],'FontSize',12);
for j=1:sum(p_adjM<0.05)
    text(j, max(obs_filt(j),exp_filt(j))+motif_ylim/30, relative_filt_str(j), 'HorizontalAlignment', 'center','FontSize',13);
    text(j, max(obs_filt(j),exp_filt(j))+motif_ylim/15, signifM_filt(j), 'HorizontalAlignment', 'center','FontSize',13,'FontWeight','bold');
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
end

%volcano plot
if min(pre_relative_obs)==0
    minrelobs=min(pre_relative_obs(pre_relative_obs>0));
else
    minrelobs=0;
end

relative_obs_volcano=log2(pre_relative_obs+minrelobs);
p_volcano=-log10(p_adjM+min(p_adjM(p_adjM~=0)));

pre_text_x=relative_obs_volcano(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));
pre_text_y=p_volcano(pre_relative_obs>=cut_relativeobs | pre_relative_obs<=(1/cut_relativeobs));

if volcanodotlabel==1
    if volcanoplotsize==0
        figure;scatter(relative_obs_volcano,p_volcano,50,'b', 'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0);hold on;
        yline(-log10(0.05),'--r','LineWidth',1);
        xline(log2(1),'--r','LineWidth',1);
        text(pre_text_x(p_adjM<0.05),pre_text_y(p_adjM<0.05),string(1:sum(p_adjM<0.05)),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',15);
        hold off;
        if sum(isinf([log2(relative_obs) log10(p_adjM)]))>0 
            set(gca,'FontSize',13,'XLim',[min(relative_obs_volcano) max(relative_obs_volcano)],'YLim',[min(p_volcano) max(p_volcano)]);
        else
            set(gca,'FontSize',13);
        end
    else
        figure;scatter(relative_obs_volcano,p_volcano,50,'b', 'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0);hold on;
        yline(-log10(0.05+min(p_adjM)),'--r','LineWidth',1);
        xline(log2(1),'--r','LineWidth',1);
        text(pre_text_x(p_adjM<0.05),pre_text_y(p_adjM<0.05),string(1:sum(p_adjM<0.05)),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',15);
        hold off;
        if sum(isinf([log2(relative_obs) log10(p_adjM)]))>0 
            set(gca,'FontSize',13);
            xlim([volcanoplotsize(1) volcanoplotsize(2)]);
            ylim([volcanoplotsize(3) volcanoplotsize(4)]);%volcanoplotsize determine x and y limits
        else
            set(gca,'FontSize',13);
            xlim([volcanoplotsize(1) volcanoplotsize(2)]);
            ylim([volcanoplotsize(3) volcanoplotsize(4)]);
        end
    end
elseif volcanodotlabel==0
    if volcanoplotsize==0
        figure;scatter(relative_obs_volcano,p_volcano,50,'b', 'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0);hold on;
        yline(-log10(0.05+min(p_adjM)),'--r','LineWidth',1);
        xline(log2(1),'--r','LineWidth',1);
        hold off;
        set(gca,'FontSize',13);
    else
        figure;scatter(relative_obs_volcano,p_volcano,50,'b', 'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0);hold on;
        yline(-log10(0.05+min(p_adjM)),'--r','LineWidth',1);
        xline(log2(1),'--r','LineWidth',1);
        hold off;
        set(gca,'FontSize',13);
        xlim([volcanoplotsize(1) volcanoplotsize(2)]);
        ylim([volcanoplotsize(3) volcanoplotsize(4)]);
    end
end


%%individual neurons of each motif plot (normalized)
mark_mtf="";
for i=1:size(mark_targets_filt,2)
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
    writematrix(Pconnected_motif_aslabel{i},'coprojection and motif_'+string(cutoff)+'.xlsx','Sheet','Sheet1','Range',("B"+num2str(i)));
end

save('motif_'+string(cutoff)+'.mat','matrix_binary','obs','exp','p_adjM','signifM','relative_obs','Pconnected_motif_asidx','Pconnected_motif_aslabel','norm_neurons_mtf_filt','norm_neurons_mtfonly_filt');
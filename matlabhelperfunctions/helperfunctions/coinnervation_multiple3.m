function coinner_multiple=coinnervation_multiple3(matrix,region_names,cutoff,colorscale,dicescale,min_region,min_cts,pvaluecut,N_t_input)

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 


%binarize projection matrix
matrix_binary_tmp=matrix>cutoff;
matrix_binary_tmp2=matrix_binary_tmp(:,sum(matrix_binary_tmp,1)>=min_region);%enforce a minimum number of cells projecting to each analyses area
matrix_binary=matrix_binary_tmp2(sum(matrix_binary_tmp2,2)~=0,:);
region_names=region_names(sum(matrix_binary_tmp,1)>=min_region);


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


nTargets = size(matrix_binary,2);

    % 각 target의 hit 확률
    p = sum(matrix_binary, 1) / N_t;

    coinner_multiple = cell(nTargets, 1);

    for k = 2:nTargets

        combs = nchoosek(1:nTargets, k);
        nComb = size(combs, 1);

        % region name 조합 생성
        RegionNames = strings(nComb, 1);
        for i = 1:nComb
            RegionNames(i) = strjoin(region_names(combs(i,:)), "-");
        end

        % joint probability
        p_comb = prod(p(combs), 2);

        % 기대값
        Exp = N_t * p_comb;

        % 관찰값
        Obs = zeros(nComb,1);
        for i = 1:nComb
            cols = combs(i,:);
            Obs(i) = sum(all(matrix_binary(:, cols) == 1, 2));
        end

        % over/under ratio
        overunder = Obs ./ Exp;

        % binomial test p-value
        pval = zeros(nComb,1);
        for i = 1:nComb
            if Obs(i)<=Exp(i)
                pval(i) = binocdf(Obs(i)-1, N_t, p_comb(i));
            else
                pval(i) = 1 - binocdf(Obs(i)-1, N_t, p_comb(i));
            end
        end


        % Bonferroni correction
        pval_bonf = min(pval * nComb, 1);

        % significance labels
        significance = strings(nComb,1);
        for i = 1:nComb
            if pval_bonf(i) < 1e-4
                significance(i) = "****";
            elseif pval_bonf(i) < 1e-3
                significance(i) = "***";
            elseif pval_bonf(i) < 1e-2
                significance(i) = "**";
            elseif pval_bonf(i) < 5e-2
                significance(i) = "*";
            else
                significance(i) = "NS";
            end
        end

        % 최종 테이블 (Obs → overunder → pval 순서 유지)
        T = table(...
            combs, RegionNames, p_comb, Exp, Obs, overunder, pval, pval_bonf, significance, ...
            'VariableNames', ...
            {'RegionIdx','RegionNames','P_comb','Exp','Obs','overunder','pval','pval_bonf','significance'});

        coinner_multiple{k} = T;
    end

    marktargets={};
    for i=2:nTargets
        regionidx_tmp=coinner_multiple{i}.RegionIdx;
        marktargets{i}=zeros(nTargets,size(regionidx_tmp,1));
        for j=1:size(regionidx_tmp,1)
             marktargets{i}(regionidx_tmp(j,:),j)=1;
        end
    end

%%plot 
% significant only
for i=2:nTargets
    
    pre_obs=coinner_multiple{i}.Obs';
    pre_exp=coinner_multiple{i}.Exp';
    pre_relative=coinner_multiple{i}.overunder';
    pre_p_adjM=coinner_multiple{i}.pval_bonf';
    pre_signifM=coinner_multiple{i}.significance;
    pre_plot_target=pre_relative.*marktargets{i};
    pre_plot_target(marktargets{i}==0)=1;
    label_target_all=repmat(region_names',1,size(pre_plot_target,2));
    pre_label_target=strings(size(pre_plot_target));
    pre_label_target(logical(marktargets{i}))=label_target_all(logical(marktargets{i}));
    
    %filt by min_cts
    obs=pre_obs(pre_obs>min_cts | pre_exp>min_cts);
    exp=pre_exp(pre_obs>min_cts | pre_exp>min_cts);
    relative=pre_relative(pre_obs>min_cts | pre_exp>min_cts);
    p_adjM=pre_p_adjM(pre_obs>min_cts | pre_exp>min_cts);
    signifM=pre_signifM(pre_obs>min_cts | pre_exp>min_cts);
    plot_target=pre_plot_target(:,pre_obs>min_cts | pre_exp>min_cts);
    label_target=pre_label_target(:,pre_obs>min_cts | pre_exp>min_cts);
   
    %filt by significance
     
    obs_filt=obs(p_adjM<pvaluecut);
    exp_filt=exp(p_adjM<pvaluecut);
    relative_filt=relative(p_adjM<pvaluecut);
    relative_filt_str=string(round(relative_filt,2));
    relative_filt_str(relative_filt>100)=">100";
    signifM_filt=signifM(p_adjM<pvaluecut);
    plot_target_filt=plot_target(:,p_adjM<pvaluecut);
    label_target_filt=label_target(:,p_adjM<pvaluecut);

    if isempty(obs_filt)==1
        continue;
    end
    
    figure;set(gcf,'Units', 'normalized', 'Position', [0, 0,sum(p_adjM<pvaluecut)*0.03, 0.5]);
    ax1=subplot(30,1,1:16);
    bar(1:length(obs_filt),[obs_filt;exp_filt]);
    motif_ylim=max([obs_filt exp_filt])*1.1;
    set(gca,'XTick',[],'YLim',[0 motif_ylim],'FontSize',12);
    for j=1:length(obs_filt)
        text(j, max(obs_filt(j),exp_filt(j))+motif_ylim/30, relative_filt_str(j), 'HorizontalAlignment', 'center','FontSize',13);
        text(j, max(obs_filt(j),exp_filt(j))+motif_ylim/15, signifM_filt(j), 'HorizontalAlignment', 'center','FontSize',13,'FontWeight','bold');
    end
    ax2=subplot(30,1,17:25);
    imagesc(plot_target_filt);clim(colorscale);colormap(redblue);set(gca,'ColorScale','log');
    for k=1:nTargets
        for j=1:sum(p_adjM<pvaluecut)
            text(j, 1+nTargets, string(j), 'HorizontalAlignment', 'center','FontSize',15);
            text(j, k, label_target_filt(k,j), 'HorizontalAlignment', 'center','FontSize',11,'FontWeight','bold');
        end
    end
    linkaxes([ax1 ax2],'x');
    set(gca,'XTick',[],'YTick',[]);
    filename = sprintf('%d_targets'+string(cutoff)+'.pdf', i);  % 01,02처럼 두 자리로 저장
    exportgraphics(gcf, filename, 'ContentType', 'vector');
    exportgraphics(gcf, '%d_targets'+string(cutoff)+'.png');
    close all;

figure;set(gcf,"Position",[0 0 60*length(Plabels) 40*length(Plabels)]);

ax = gca;ax.YTick=[];ax.XTick=[];
axis off;
hold on;
for i = 1:length(conditionalp_dice)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

    figure;set(gcf,'Units', 'normalized', 'Position', [0, 0,sum(p_adjM<pvaluecut)*0.03, 0.5]);
    ax1=subplot(30,1,1:16);
    imagesc(obs);clim([0 dicescale]);
    for k=1:nTargets
        for j=1:sum(p_adjM<pvaluecut)
            text(j, 1+nTargets, string(j), 'HorizontalAlignment', 'center','FontSize',15);
            text(j, k, label_target_filt(k,j), 'HorizontalAlignment', 'center','FontSize',11,'FontWeight','bold');
        end
    end
    linkaxes([ax1 ax2],'x');
    set(gca,'XTick',[],'YTick',[]);
    filename = sprintf('%d_targets'+string(cutoff)+'.pdf', i);  % 01,02처럼 두 자리로 저장
    exportgraphics(gcf, filename, 'ContentType', 'vector');
    exportgraphics(gcf, '%d_targets'+string(cutoff)+'.png');
    close all;
end

save("coinner_multiple"+string(cutoff)+".mat",'coinner_multiple');
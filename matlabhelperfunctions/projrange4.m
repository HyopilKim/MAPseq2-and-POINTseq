function projrange4(data1,data2,sourceColors,targets)
data1=data1(sum(data1,2)>0,:);
data2=data2(sum(data2,2)>0,:);

aligned1_sorted=zeros(size(data1,1),length(targets));
aligned2_sorted=zeros(size(data2,1),length(targets));
[~,idx1]=max(data1,[],2);
[~,idx2]=max(data2,[],2);

aligned1_cell={};
aligned2_cell={};

for i=1:size(data1,1)
    aligned1_sorted(i,:)=sort(data1(i,:),'descend');
    for j=1:size(data1,2)
        aligned1_cell{j}=data1(idx1==j,:);
    end
    aligned1_cell{size(data1,2)+1}=aligned1_sorted;
end

for i=1:size(data2,1)
    aligned2_sorted(i,:)=sort(data2(i,:),'descend');
    for j=1:size(data2,2)
        aligned2_cell{j}=data2(idx2==j,:);
    end
    aligned2_cell{size(data2,2)+1}=aligned2_sorted;
end

% mean and SEM
upper1={};
lower1={};
mean1={};
mean2={};
upper2={};
lower2={};
sem1={};
sem2={};
for i=1:length(aligned1_cell)
    mean1{i}=mean(aligned1_cell{i}, 1); 
    sem1{i} =std(aligned1_cell{i}, 0, 1) / sqrt(size(aligned1_cell{i}, 1));
    upper1{i} = mean1{i} + sem1{i}; 
    lower1{i} = mean1{i} - sem1{i}; 
end

for i=1:length(aligned2_cell)
    mean2{i}=mean(aligned2_cell{i}, 1); 
    sem2{i} =std(aligned2_cell{i}, 0, 1) / sqrt(size(aligned2_cell{i}, 1));
    upper2{i} = mean2{i} + sem2{i}; 
    lower2{i} = mean2{i} - sem2{i}; 
end

rm_p={};
rm_p_inter=[];
bonferroni=ones(length(targets)+1,length(targets));
region_idx=zeros(length(targets)+1,length(targets));

for i=1:length(aligned1_cell)
    if ~isempty(aligned1_cell{i}) && ~isempty(aligned2_cell{i})
        pre_group1=aligned1_cell{i};pre_group2=aligned2_cell{i};
        num_regions=sum((sum(pre_group1,1)>0 & sum(pre_group2,1)>0));% compare where not all the values are zero in a group
        group1=pre_group1(:,(sum(pre_group1,1)>0 & sum(pre_group2,1)>0));
        group2=pre_group2(:,(sum(pre_group1,1)>0 & sum(pre_group2,1)>0));
        region_idx(i,:)=(sum(pre_group1,1)>0 & sum(pre_group2,1)>0);
        anova_data = [group1; group2]; 
        group_factor = [ones(size(group1, 1), 1); 2 * ones(size(group2, 1), 1)];
        anova_data = [group_factor anova_data]; 
        region_factor = repmat(1:num_regions, size(anova_data, 1), 1);
        region_factor = region_factor(:);
        Regions=cellstr(targets((sum(pre_group1,1)>0 & sum(pre_group2,1)>0)));
        T = array2table(anova_data,VariableNames=['Source', Regions]);
        T.Source=string(T.Source);
        rm=fitrm(T, strcat(Regions{1},'-',Regions{end},' ~ Source'));
        ranova_results = ranova(rm);
        rm_p{i}=ranova_results;
        rm_p_inter(i)=ranova_results.pValueGG(2);
        posthoc_results = table2array(multcompare(rm,'Source','By','Time', 'ComparisonType', 'bonferroni'));
        posthoc_results = str2double(posthoc_results(1:2:end,6));
        bonferroni(i,sum(pre_group1,1)>0 & sum(pre_group2,1)>0)=posthoc_results;
    end
end

signif=repmat("",length(targets)+1,length(targets));
signif(bonferroni>=0.01 & bonferroni<0.05)="*";
signif(bonferroni>=0.001 & bonferroni<0.01)="**";
signif(bonferroni>=0.0001 & bonferroni<0.001)="***";
signif(bonferroni<0.0001)="****";

%%plot
figure; set(gcf,'Units', 'normalized', 'Position', [0, 0, 0.2*length(targets)/3, 0.6]);findfigs;
for i=1:length(aligned1_cell)  
    clf;
    hold on; 
    fill([1:length(targets), fliplr(1:length(targets))], [upper1{i}, fliplr(lower1{i})], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(1:length(targets), mean1{i}, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
    fill([1:length(targets), fliplr(1:length(targets))], [upper2{i}, fliplr(lower2{i})], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(1:length(targets), mean2{i}, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
    for j=1:length(targets)
        if rm_p_inter(i)<0.05
            text(j,max(mean1{i}(j),mean2{i}(j))+0.05,signif(i,j),"FontSize",28,"FontWeight",'bold','HorizontalAlignment','center')
        end
    end
    hold off;
    if i<length(aligned1_cell) 
        xticks(1:length(targets));set(gca,'FontSize',40,'YLim',[0 1],'XTickLabel',targets);
        pause(0.5);exportgraphics(gcf,(targets(i)+"_projdistribution.png"));
    else
        xticks(1:length(targets));set(gca,'FontSize',40,'YLim',[0 1],'XTickLabel',string(1:length(targets)));
        pause(0.5);exportgraphics(gcf,"all_projdistribution.png");
    end
end


mean1_grid={};
for i=1:size(data1,2)
    tmp=NaN(4,4);
    means=mean1{i};
    tmp(1,1:3)=means(1:3);
    tmp(2,1:3)=means(4:6);
    tmp(3,1:2)=means(7:8);
    tmp(4,4)=means(9);
    mean1_grid{i}=tmp;
end

mean2_grid={};
for i=1:size(data2,2)
    tmp=NaN(4,4);
    means=mean2{i};
    tmp(1,1:3)=means(1:3);
    tmp(2,1:3)=means(4:6);
    tmp(3,1:2)=means(7:8);
    tmp(4,4)=means(9);
    mean2_grid{i}=tmp;
end



%figure;set(gcf,'Units', 'normalized', 'Position', [0 0 size(mean1_grid{i},2)*0.1, size(mean1_grid{i},1)*0.1]); findfigs;
%for i=1:size(data1,2)  
%    clf;
%    hold on; 
%    h=imagesc(mean1_grid{i}(end:-1:1,:));colormap("parula"); 
%    hold off;
%    set(gca,'FontSize',19,'XTick',[],'YTick',[]);
%    alphaData = ~isnan(mean1_grid{i}(end:-1:1,:));
%    set(h, 'AlphaData', alphaData);
%    pause(1);exportgraphics(gcf,(targets(i)+"_projdistribution1.png"));
%end

%figure;set(gcf,'Units', 'normalized', 'Position', [0 0 size(mean2_grid{i},2)*0.1, size(mean2_grid{i},1)*0.1]); findfigs;
%for i=1:size(data2,2)  
%    clf;
%    hold on; 
%    h=imagesc(mean2_grid{i}(end:-1:1,:));colormap("parula"); 
%    hold off;
%    set(gca,'FontSize',19,'XTick',[],'XTickLabelRotation',45,'YTick',[]);
%    alphaData = ~isnan(mean2_grid{i}(end:-1:1,:));
%    set(h, 'AlphaData', alphaData);
%    pause(1);exportgraphics(gcf,(targets(i)+"_projdistribution2.png"));
%end






save("projection_range.mat","aligned1_cell","aligned2_cell","mean1_grid","mean2_grid","bonferroni","rm_p","rm_p_inter");
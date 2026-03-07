function aligned_cell=projrange5(data,targets)

data=data(sum(data,2)>0,:);
data_rd=data(randperm(size(data,1)),:);
[~,maxidx]=max(data_rd,[],2);

aligned=zeros(size(data_rd));
aligned_cell={};
mean_aligned=zeros(size(data_rd,2),size(data_rd,2));
k=1;
for i=1:length(targets)
    num_i_neurons=sum(maxidx==i);
    aligned(k:k+num_i_neurons-1,:)=data_rd((maxidx==i),:);
    mean_aligned(i,:)=mean(data_rd((maxidx==i),:),1);
    aligned_cell{i}=data_rd((maxidx==i),:);
    k=k+num_i_neurons;
end

% plotting 
figure;set(gcf,'Units', 'normalized', 'Position', [0, 0, 0.4, 0.8]);
imagesc(aligned);colormap('parula');
set(gca,'xtick',1:numel(targets),'YAxisLocation', 'right','XTickLabel',targets,'XTickLabelRotation',90)
ax = gca;ax.YAxis.FontSize = 16;ax.XAxis.FontSize = 18;

figure;set(gcf,'Units', 'normalized', 'Position', [0, 0, 0.4, 0.3]);
imagesc(mean_aligned);colormap('parula');
set(gca,'xtick',1:numel(targets),'YAxisLocation', 'right','XTickLabel',targets,'XTickLabelRotation',90)
ax = gca;ax.YAxis.FontSize = 16;ax.XAxis.FontSize = 18;

findfigs;
save("projection_maxaligned.mat","aligned","mean_aligned","maxidx","aligned_cell");
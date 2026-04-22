function projrange3(data1,data2,sourceColors,num_range)
data1=data1(sum(data1,2)>0,:);
data2=data2(sum(data2,2)>0,:);
data1=data1./repmat(sum(data1,2),1,size(data1,2));
data2=data2./repmat(sum(data2,2),1,size(data2,2));

aligned_r1=zeros(size(data1,1),num_range);
aligned_r2=zeros(size(data2,1),num_range);
aligned_c1=zeros(size(data1,1),num_range);
aligned_c2=zeros(size(data2,1),num_range);
[~,idx1]=max(data1,[],2);
[~,idx2]=max(data2,[],2);

for i=1:size(data1,1)
    if idx1(i)<=(size(data1,2)-(num_range-1))
        for j=1:num_range
            aligned_r1(i,j)=data1(i,idx1(i)+j-1);
        end
    end
    if idx1(i)>=num_range
        for j=1:num_range
            aligned_c1(i,j)=data1(i,(idx1(i)-num_range+j));
        end
    end
    aligned_r1=aligned_r1(max(aligned_r1,[],2)>0,:);
    aligned_c1=aligned_c1(max(aligned_c1,[],2)>0,:);
end

for i=1:size(data2,1)
    if idx2(i)<=(size(data2,2)-(num_range-1))
        for j=1:num_range
            aligned_r2(i,j)=data2(i,idx2(i)+j-1);
        end
    end
    if idx2(i)>=num_range
        for j=1:num_range
            aligned_c2(i,j)=data2(i,(idx2(i)-num_range+j));
        end
    end
    aligned_r2=aligned_r2(max(aligned_r2,[],2)>0,:);
    aligned_c2=aligned_c2(max(aligned_c2,[],2)>0,:);
end

% mean and SEM
mean_r1= mean(aligned_r1, 1); 
sem_r1 = std(aligned_r1, 0, 1) / sqrt(size(aligned_r1, 1)); 
mean_r2= mean(aligned_r2, 1); 
sem_r2 = std(aligned_r2, 0, 1) / sqrt(size(aligned_r2, 1)); 
mean_c1= mean(aligned_c1, 1); 
sem_c1 = std(aligned_c1, 0, 1) / sqrt(size(aligned_c1, 1)); 
mean_c2= mean(aligned_c2, 1); 
sem_c2 = std(aligned_c2, 0, 1) / sqrt(size(aligned_c2, 1)); 


% upper and lower limits by SEM
upper_r1 = mean_r1 + sem_r1; 
lower_r1 = mean_r1 - sem_r1; 
upper_r2 = mean_r2 + sem_r2; 
lower_r2 = mean_r2 - sem_r2; 
upper_c1 = mean_c1 + sem_c1; 
lower_c1 = mean_c1 - sem_c1; 
upper_c2 = mean_c2 + sem_c2; 
lower_c2 = mean_c2 - sem_c2; 


% X range 
x_c = zeros(1,num_range);
for i=1:(num_range-1)
    x_c(num_range-i)=(-0.3*i);
end

x_r = zeros(1,num_range);
for i=1:num_range
    x_r(i)=0.3*(i-1);
end


%% r plot
figure; set(gcf,'Units', 'normalized', 'Position', [0, 0, 0.2*num_range/3, 0.6]);hold on;
fill([x_r, fliplr(x_r)], [upper_r1, fliplr(lower_r1)], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_r, mean_r1, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
fill([x_r, fliplr(x_r)], [upper_r2, fliplr(lower_r2)], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_r, mean_r2, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
findfigs;
xticks(x_r);set(gca,'XTickLabel',string(x_r),'FontSize',18,'YLim',[0 1]);




%% c plot
figure; set(gcf,'Units', 'normalized', 'Position', [0, 0, 0.2*num_range/3, 0.6]);hold on;
fill([x_c, fliplr(x_c)], [upper_c1, fliplr(lower_c1)], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_c, mean_c1, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
fill([x_c, fliplr(x_c)], [upper_c2, fliplr(lower_c2)], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_c, mean_c2, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
findfigs;
xticks(x_c);set(gca,'XTickLabel',string(x_c),'FontSize',18,'YLim',[0 1]);


save("projection_range.mat","aligned_r1","aligned_c1","aligned_r2","aligned_c2");
function [aligned_r1,aligned_c1,aligned_r2,aligned_c2,aligned1,aligned2]=projrange(data1,data2,sourceColors)
data1=data1(sum(data1,2)>0,:);
data2=data2(sum(data2,2)>0,:);
data1=data1./repmat(sum(data1,2),1,size(data1,2));
data2=data2./repmat(sum(data2,2),1,size(data2,2));

forrange=min(length(data1),lenth(data2));

if mod(forrange,2)==0
    proj_range=forrange-1;
else
    proj_range=

aligned1=zeros(size(data1,1),7);
aligned2=zeros(size(data2,1),7);
aligned_r1=zeros(size(data1,1),5);
aligned_r2=zeros(size(data2,1),5);
aligned_c1=zeros(size(data1,1),5);
aligned_c2=zeros(size(data2,1),5);
[~,idx1]=max(data1,[],2);
[~,idx2]=max(data2,[],2);

for i=1:size(data1,1)
    if (idx1(i)>=4) && (idx1(i)<=(size(data1,2)-3))
        aligned1(i,7)=data1(i,idx1(i)+3);
        aligned1(i,6)=data1(i,idx1(i)+2);
        aligned1(i,5)=data1(i,idx1(i)+1);
        aligned1(i,4)=data1(i,idx1(i));
        aligned1(i,3)=data1(i,idx1(i)-1);
        aligned1(i,2)=data1(i,idx1(i)-2);
        aligned1(i,1)=data1(i,idx1(i)-3);
    end
    if idx1(i)>=5
        aligned_r1(i,5)=data1(i,idx1(i));
        aligned_r1(i,4)=data1(i,idx1(i)-1);
        aligned_r1(i,3)=data1(i,idx1(i)-2);
        aligned_r1(i,2)=data1(i,idx1(i)-3);
        aligned_r1(i,1)=data1(i,idx1(i)-4);
    end
    if idx1(i)<=(size(data1,2)-4)
        aligned_c1(i,5)=data1(i,idx1(i)+4);
        aligned_c1(i,4)=data1(i,idx1(i)+3);
        aligned_c1(i,3)=data1(i,idx1(i)+2);
        aligned_c1(i,2)=data1(i,idx1(i)+1);
        aligned_c1(i,1)=data1(i,idx1(i));
    end
    aligned1=aligned1(max(aligned1,[],2)>0,:);
    aligned_r1=aligned_r1(max(aligned_r1,[],2)>0,:);
    %aligned_r1=aligned_r1./repmat(sum(aligned_r1,2),1,size(aligned_r1,2));
    aligned_c1=aligned_c1(max(aligned_c1,[],2)>0,:);
    %aligned_c1=aligned_c1./repmat(sum(aligned_c1,2),1,size(aligned_c1,2));
end

for i=1:size(data2,1)
    if (idx2(i)>=4) && (idx2(i)<=(size(data2,2)-3))
        aligned2(i,7)=data2(i,idx2(i)+3);
        aligned2(i,6)=data2(i,idx2(i)+2);
        aligned2(i,5)=data2(i,idx2(i)+1);
        aligned2(i,4)=data2(i,idx2(i));
        aligned2(i,3)=data2(i,idx2(i)-1);
        aligned2(i,2)=data2(i,idx2(i)-2);
        aligned2(i,1)=data2(i,idx2(i)-3);
    end
    if idx2(i)>=5
        aligned_r2(i,5)=data2(i,idx2(i));
        aligned_r2(i,4)=data2(i,idx2(i)-1);
        aligned_r2(i,3)=data2(i,idx2(i)-2);
        aligned_r2(i,2)=data2(i,idx2(i)-3);
        aligned_r2(i,1)=data2(i,idx2(i)-4);
    end
    if idx2(i)<=(size(data2,2)-4)
        aligned_c2(i,5)=data2(i,idx2(i)+4);
        aligned_c2(i,4)=data2(i,idx2(i)+3);
        aligned_c2(i,3)=data2(i,idx2(i)+2);
        aligned_c2(i,2)=data2(i,idx2(i)+1);
        aligned_c2(i,1)=data2(i,idx2(i));
    end
    aligned2=aligned2(max(aligned2,[],2)>0,:);
    aligned_r2=aligned_r2(max(aligned_r2,[],2)>0,:);
    %aligned_r2=aligned_r2./repmat(sum(aligned_r2,2),1,size(aligned_r2,2));
    aligned_c2=aligned_c2(max(aligned_c2,[],2)>0,:);
    %aligned_c2=aligned_c2./repmat(sum(aligned_c2,2),1,size(aligned_c2,2));
end

% 평균 및 SEM 계산
mean1= mean(aligned1, 1); 
sem1 = std(aligned1, 0, 1) / sqrt(size(aligned_r1, 1)); 
mean2= mean(aligned2, 1); 
sem2 = std(aligned2, 0, 1) / sqrt(size(aligned_r2, 1)); 
mean_r1= mean(aligned_r1, 1); 
sem_r1 = std(aligned_r1, 0, 1) / sqrt(size(aligned_r1, 1)); 
mean_r2= mean(aligned_r2, 1); 
sem_r2 = std(aligned_r2, 0, 1) / sqrt(size(aligned_r2, 1)); 
mean_c1= mean(aligned_c1, 1); 
sem_c1 = std(aligned_c1, 0, 1) / sqrt(size(aligned_c1, 1)); 
mean_c2= mean(aligned_c2, 1); 
sem_c2 = std(aligned_c2, 0, 1) / sqrt(size(aligned_c2, 1)); 


% 상한/하한 계산
upper1 = mean1 + sem1; % 평균 + SEM
lower1 = mean1 - sem1; % 평균 - SEM
upper2 = mean2 + sem2; % 평균 + SEM
lower2 = mean2 - sem2; % 평균 - SEM
upper_r1 = mean_r1 + sem_r1; % 평균 + SEM
lower_r1 = mean_r1 - sem_r1; % 평균 - SEM
upper_r2 = mean_r2 + sem_r2; % 평균 + SEM
lower_r2 = mean_r2 - sem_r2; % 평균 - SEM
upper_c1 = mean_c1 + sem_c1; % 평균 + SEM
lower_c1 = mean_c1 - sem_c1; % 평균 - SEM
upper_c2 = mean_c2 + sem_c2; % 평균 + SEM
lower_c2 = mean_c2 - sem_c2; % 평균 - SEM


% X축 값 
x=[-0.9,-0.6, -0.3, 0, 0.3, 0.6, 0.9];
x_r = [-1.2, -0.9, -0.6, -0.3, 0];
x_c = [0, 0.3, 0.6, 0.9, 1.2];


%% r+c plot
figure; hold on;
fill([x, fliplr(x)], [upper1, fliplr(lower1)], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, mean1, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
fill([x, fliplr(x)], [upper2, fliplr(lower2)], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, mean2, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
findfigs;
% X축 설정
xticks(x);
xticklabels({'-0.9','-0.6', '-0.3', '0','0.3','0.6','0.9'});

xlabel('Relative position');
ylabel('Projection strength');


%% r 플롯
figure; hold on;
fill([x_r, fliplr(x_r)], [upper_r1, fliplr(lower_r1)], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_r, mean_r1, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
fill([x_r, fliplr(x_r)], [upper_r2, fliplr(lower_r2)], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_r, mean_r2, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
findfigs;
% X축 설정
xticks(x_r);
xticklabels({'-1.2','-0.9','-0.6', '-0.3', '0'});

xlabel('Relative position');
ylabel('Projection strength');



%% c 플롯
figure; hold on;
fill([x_c, fliplr(x_c)], [upper_c1, fliplr(lower_c1)], [1 0.5 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_c, mean_c1, '-o','Color',sourceColors(1,:), 'LineWidth', 1.5);
fill([x_c, fliplr(x_c)], [upper_c2, fliplr(lower_c2)], [0.5 0.6 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_c, mean_c2, '-o','Color',sourceColors(2,:), 'LineWidth', 1.5);
findfigs;
% X축 설정
xticks(x_c);
xticklabels({'0', '0.3', '0.6', '0.9', '1.2'});

xlabel('Relative position');
ylabel('Projection strength');



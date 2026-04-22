function quantproj(matrix,nperm,slim,dscale,dolabel)
%% shuffled
D = zeros(size(matrix,2));
for i = 1:size(matrix,2)
    for j =1:size(matrix,2)
        d=[];
        for k=1:size(matrix,1)
            if (matrix(k, i) + matrix(k, j))>0
                d=[d,abs(matrix(k, i) - matrix(k, j))/(matrix(k, i) + matrix(k, j))];
            end
        end
        D(i,j)=mean(d);
    end
end
similarity=1-D;


similarity_rd=zeros(size(similarity,1),size(similarity,2),nperm);
for n=1:nperm
    rd=shuffincol(matrix);
    D = zeros(size(rd,2));
    for i = 1:size(rd,2)
        for j =1:size(rd,2)
            d=[];
            for k=1:size(rd,1)
                if (rd(k, i) + rd(k, j))>0
                    d=[d,abs(rd(k, i) - rd(k, j))/(rd(k, i) + rd(k, j))];
                end
            end
            D(i,j)=mean(d);
        end
    end
    similarity_rd(:,:,n)=1-D;
end
mean_similarity_rd=mean(similarity_rd,3);


figure;set(gcf,"Position",[0 0 60*size(matrix,2) 40*size(matrix,2)]);
imagesc(similarity);clim([0 slim])
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(D)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

figure;set(gcf,"Position",[0 0 60*size(matrix,2) 40*size(matrix,2)]);
imagesc(mean_similarity_rd);clim([0 slim])
ax = gca;
ax.YTick=[];
ax.XTick=[];
axis off;
hold on;
for i = 1:length(D)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;



[X, Y, N] = size(similarity_rd);
p_values = nan(X, Y);

for i = 1:X
    for j = 1:Y
        null_vals = squeeze(similarity_rd(i, j, :));
        obs_val = similarity(i, j);
        null_mean = mean(null_vals);
        p = (sum(abs(null_vals - null_mean) >= abs(obs_val - null_mean)) + 1) / (N + 1);
        p_values(i, j) = p;
    end
end



% Apply Bonferroni correction
p_adjP = p_values * nchoosek(size(similarity,2),2);
signifP=strings(size(p_adjP));
signifP(p_adjP>=0.05)="NS";
signifP(p_adjP<0.05 & p_adjP>=0.01)="*";
signifP(p_adjP<0.01 & p_adjP>=0.001)="**";
signifP(p_adjP<0.001 & p_adjP>=0.0001)="***";
signifP(p_adjP<0.0001)="****";
signifP(p_adjP==1)="";

depen=similarity./mean_similarity_rd;

n = 256;  % number of colors
redblue = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); 
        ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)']; 

figure;imagesc(depen);clim(dscale);colormap(redblue);
set(gcf,"Position",[0 0 60*size(matrix,2) 40*size(matrix,2)]);
set(gca,'ColorScale','log');
ax = gca;
ax.YTick=[];
ax.XTick=[];
[rows, cols] = size(depen);
for i = 1:rows
   for j = 1:cols
        if dolabel==1
        text(j, i, num2str(depen(i,j),'%.2f'), 'HorizontalAlignment', 'center','FontSize',10);
        text(j, i-0.3, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',10);
        else
        text(j, i, signifP(i,j), 'HorizontalAlignment', 'center','FontSize',13);
        end
    end
end
hold on;
for i = 1:length(depen)
    rectangle('Position', [i-0.5, i-0.5, 1, 1], 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
end
hold off;

findfigs;

save('quantprojoutput.mat','depen','similarity_rd');
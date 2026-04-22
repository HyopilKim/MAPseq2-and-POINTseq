function [percent, total]=numproj_injthresh(barcodematrix,xlimit,projthresh,sourcesite,projsite,n)

% barcodematrix: input raw matrix
% xlimit is x axis limit to display
% projthresh is threshold for target sites UMI
% sourcesite is the number assigned to source regions
% projsite is the number assigned to target regions
% n is maximum number of targets to display

percent=zeros(xlimit,n);
total=1:xlimit;
for i=1:xlimit
    sourcethresh=i;
    B=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
    B_tar=B(:,projsite);
    B_tar(B_tar>0)=1;
    num=sum(B_tar,2);
    num(num>n)=n;
    binEdges = 1:(n+1);
    count=histcounts(num,binEdges);
    total(i)=sum(count);
    percent(i,:)=count/total(i);
end

b=bar(percent,'stacked','FaceColor','flat');
colors=parula(n);
for k = 1:n
    b(k).FaceColor = 'flat';
    b(k).EdgeColor = 'none';
    b(k).CData = colors(k,:);  
    b(k).BarWidth=1;
end
ylim([0 1]);

figure;scatter(1:xlimit,total);
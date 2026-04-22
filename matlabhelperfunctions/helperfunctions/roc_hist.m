function [true_cumul,true_bin, num]=roc_hist(sourceratio,ratio,sourceidx,binwidth)

edges = 0:binwidth:1;



ratio_only=ratio(:,1);
idx=ratio(:,2);
true_cumul=zeros(length(1:(1/binwidth)),1);
x=zeros(length(1:(1/binwidth)),1);

for i=1:(1/binwidth)
    r=binwidth*i;
    true=sum(idx(ratio_only<=r)==sourceidx)/sum(ratio_only<=r);
    true_cumul(i)=true;
    x(i)=r;
end
roc_plot=[x,true_cumul];



% create bins

[~,~,idx]=histcounts(ratio(:,1),edges);


%in each bin, calculate accuracy
true_bin=zeros(length(edges)-1,1);
num=zeros(length(edges)-1,1);
for i=1:(length(edges)-1)
    r=ratio(idx==i,1);
    s=ratio(idx==i,2);
    true_bin(i)=sum((s==sourceidx))/length(r);
    num(i)=length((s==sourceidx));
end


binCenters = edges(1:end-1) + diff(edges)/2;
figure;plot(edges(2:end),roc_plot(:,2),'-o','LineWidth',1.5);ylim([0 1]);
figure;bar(binCenters,histcounts(sourceratio,edges)/sum(histcounts(sourceratio,edges)),'BarWidth',1);ylim([0 1]);
figure;bar(binCenters,true_bin,'BarWidth',1);ylim([0 1]);
figure;bar(binCenters,num,'BarWidth',1,'FaceColor','#FFA500');
findfigs;

function [true_cumul,true_bin, num]=roc_hist(sourceratio,ratio,binwidth)

%draw cumulative accuracy, accuracy in each bin, number of neurons in each bin
%ratio: maxsource vs maxtarget ratio



edges = 0:binwidth:1;



ratio_only=ratio(:,1);
idx=ratio(:,2);
true_cumul=zeros(length(1:(1/binwidth)),1);
x=zeros(length(1:(1/binwidth)),1);

for i=1:(1/binwidth)
    r=binwidth*i;
    true=sum(idx(ratio_only<=r)==1)/sum(ratio_only<=r);
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
    true_bin(i)=sum((s==1))/length(r);
    num(i)=length((s==1));
end


binCenters = edges(1:end-1) + diff(edges)/2;
figure;bar(binCenters,roc_plot(:,2),'BarWidth',1);%cumulative accuracy
figure;bar(binCenters,histcounts(sourceratio,edges)/sum(histcounts(sourceratio,edges)),'BarWidth',1);ylim([0 1]);%sourceratio distribution
figure;bar(binCenters,true_bin,'BarWidth',1);ylim([0 1]);%accuracy in each bin
figure;bar(binCenters,num,'BarWidth',1,'FaceColor','#FFA500');%number of neurons in each bin
findfigs;

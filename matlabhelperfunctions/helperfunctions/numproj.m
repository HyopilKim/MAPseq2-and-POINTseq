function numproj(B,n,fontsize)

B(B>0)=1;
num=sum(B,2);
num(num>n)=n;
count=histcounts(num);
total=sum(count);
percent=count/total;
p=pie(percent);
textObjs = findobj(p, 'Type', 'text');
set(textObjs, 'FontSize', fontsize); 


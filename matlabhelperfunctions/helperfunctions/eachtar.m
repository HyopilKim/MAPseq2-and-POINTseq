function eachtar(data1,data2)

n=size(data1,2);
nofilt_1=cell(1,n);
nofilt_2=cell(1,n);
filt0_1=cell(1,n);
filt0_2=cell(1,n);

for i=1:n
    each1=data1(:,i);
    each2=data2(:,i);
    nofilt_1{i}=each1;
    nofilt_2{i}=each2;
    filt0_1{i}=each1(each1>0);
    filt0_2{i}=each2(each2>0);
    save('eachtar.mat','nofilt_1','nofilt_2','filt0_1','filt0_2');
end

 
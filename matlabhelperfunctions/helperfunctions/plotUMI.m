function plotUMI(B,rownumber,binwidth,limit)

for i=1:size(B,2)
    Bi=B(:,i);
    Bii=Bi(Bi>0);
    if size(Bii,1)<rownumber
        continue
    end
    histogram(Bii,'BinWidth',binwidth);title(int2str(i));xlim([0 limit]);findfigs;pause();
end
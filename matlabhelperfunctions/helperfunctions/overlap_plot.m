function [UMIoverlap]=overlap_plot(B1,B2,Bseq1,Bseq2,site)

for i = site
    x=B1(ismember(Bseq1,Bseq2,'rows'),i);
    y=B2(ismember(Bseq2,Bseq1,'rows'),i);
    UMIoverlap=[x,y];
    writematrix(UMIoverlap,strcat("UMIoverlap_site",string(i)));
end



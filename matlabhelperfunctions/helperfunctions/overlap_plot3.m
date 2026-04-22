function UMIoverlap=overlap_plot3(B1,B2,Bseq1,Bseq2,site)

UMIoverlap={};

for i = 1:site+1
    if i<site+1
        x=B1(ismember(Bseq1,Bseq2,'rows'),i);
        y=B2(ismember(Bseq2,Bseq1,'rows'),i);
        UMIoverlap{i}=[x y];
        writematrix(UMIoverlap{i},strcat("UMIoverlap_site",string(i)));
    else
        B1_mean=mean(B1,2);
        x=B1_mean(ismember(Bseq1,Bseq2,'rows'),1);
        B2_mean=mean(B2,2);
        y=B2_mean(ismember(Bseq2,Bseq1,'rows'),1);
        UMIoverlap{i}=[x y];
        writematrix(UMIoverlap{i},strcat("UMIoverlap_site",string(i)));
    end
end



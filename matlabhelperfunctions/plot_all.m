function [UMIoverlap]=plot_all(B1,B2,Bseq1,Bseq2,site)

for i = 1:site+1
    if i<site+1
        x=B1(ismember(Bseq1,Bseq2,'rows'),i);
        x_notiny=B1(~ismember(Bseq1,Bseq2,'rows'),i);
        y_notiny=zeros(size(x_notiny));
        y=B2(ismember(Bseq2,Bseq1,'rows'),i);
        y_notinx=B2(~ismember(Bseq2,Bseq1,'rows'),i);
        x_notinx=zeros(size(y_notinx));
        UMIoverlap=[[x;x_notiny;x_notinx] [y;y_notiny;y_notinx]];
        UMIoverlap=UMIoverlap(max(UMIoverlap,[],2)>0,:);
        writematrix(UMIoverlap,strcat("UMIoverlap_site",string(i)));
    else
        B1_mean=mean(B1,2);
        x=B1_mean(ismember(Bseq1,Bseq2,'rows'),1);
        x_notiny=B1_mean(~ismember(Bseq1,Bseq2,'rows'),1);
        y_notiny=zeros(size(x_notiny));
        B2_mean=mean(B2,2);
        y=B2_mean(ismember(Bseq2,Bseq1,'rows'),1);
        y_notinx=B2_mean(~ismember(Bseq2,Bseq1,'rows'),1);
        x_notinx=zeros(size(y_notinx));
        UMIoverlap=[[x;x_notiny;x_notinx] [y;y_notiny;y_notinx]];
        UMIoverlap=UMIoverlap(max(UMIoverlap,[],2)>0,:);
        writematrix(UMIoverlap,strcat("UMIoverlap_site",string(i)));
    end
end



function [UMIoverlap,b,R_squared,B]=overlapplot_allcomb(B1,B2,Bseq1,Bseq2)


B1_shared=B1(ismember(Bseq1,Bseq2,'rows'),:);
B2_shared=B2(ismember(Bseq2,Bseq1,'rows'),:);
B=[B1_shared B2_shared];

b=nchoosek(1:size(B,2),2);

UMIoverlap={};
R_squared=[];

for i=1:length(b)
    x=B(:,b(i,1));
    y=B(:,b(i,2));
    shared=[x,y];
    shared=shared(min(shared,[],2)>0,:);
    UMIoverlap{i}=shared;
    writematrix(UMIoverlap{i},strcat("UMIoverlap_site",string(b(i,1)),'-',string(b(i,2))));
    R = corrcoef(x, y);
    R_squared(i) = R(1, 2)^2;
end


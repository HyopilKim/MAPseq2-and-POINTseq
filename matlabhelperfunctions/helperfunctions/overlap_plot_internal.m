function [UMIoverlap]=overlap_plot_internal(B)

b=combntns(1:size(B,2),2);
UMIoverlap={};

for i = 1:length(b)
    idx1=b(i,1);
    idx2=b(i,2);
    x=B(:,idx1);
    y=B(:,idx2);
    UMIoverlap{i}=[x y];
    writematrix(UMIoverlap{i},strcat("UMIoverlap_site",string(b(i,1)),'-',string(b(i,2))));
end



%% function to compare the distances between a reference set of neurons and a testcase
function out=sampleandcomparepoints2(reference,test,repeats,distance)

out_test=zeros(repeats,1);
out_ref=zeros(repeats,1);
out_random=zeros(repeats,1);


%generate random neurons by sampling from a uniform distribution 
%log normalized as same as MAPseq1 and MAPseq2 comparison dataset)
random_tmp=rand(size(reference));
randomneurons=random_tmp./repmat(sum(random_tmp,2),1,size(random_tmp,2));
randomneurons=log(1+randomneurons*100);
randomneurons=randomneurons/max(max(randomneurons));

%calculate min distance within reference, between test and ref, between
%test and random

for i=1:repeats
    r=randi(size(reference,1));
    remains=reference([1:r-1 r+1:end],:);
    ref_r=reference(r,:);
    test_r=test(r,:);
    ref_mindist=min(pdist2(ref_r,remains,distance));
    test_mindist=min(pdist2(test_r,remains,distance));
    rand_mindist=min(pdist2(test_r,randomneurons([1:r-1 r+1:end],:),distance));
    out_ref(i)=ref_mindist;
    out_test(i)=test_mindist;
    out_random(i)=rand_mindist;
end

out=[out_ref out_test out_random];


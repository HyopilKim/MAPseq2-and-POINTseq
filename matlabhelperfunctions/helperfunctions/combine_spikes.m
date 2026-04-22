function newspikes=combine_spikes(old,new,oldspikes,newspikes)


for i=1:length(old)
    c=[];
    c2u=[];
    r=[];
    r2u=[];
    for j=old{i}
        c=[c;oldspikes(j).counts];
        c2u=[c2u;oldspikes(j).counts2u];
        r=[r;oldspikes(j).reads];
        r2u=[r2u;oldspikes(j).reads2u];
    end
    newspikes(new(i)).counts=c;
    newspikes(new(i)).counts2u=c2u;
    newspikes(new(i)).reads=r;
    newspikes(new(i)).reads2u=r2u;
end


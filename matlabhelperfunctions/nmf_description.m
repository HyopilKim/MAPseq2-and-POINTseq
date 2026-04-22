%% find k values by checking residuals or explained variance using Frobenius norm error and RSS (neuron or row based) 

%run NMF across k values
V=data;
fro=[];
normmoduleres=[];
krange=2:26;

for k=krange
    [W_pre,H_pre]=nnmf(V,k,'replicates',10);
    fro(k-1)=norm(V-W_pre*H_pre,'fro');
    RSS(k-1)=sum(sum((V-W_pre*H_pre).^2,2))/length(V); 
end


% values for normalization
totalfronorm=norm(V,'fro');
TSS=sum(sum((repmat(mean(V),size(V,1),1)-V).^2,2))/length(V);

% draw plots
figure;
plot(krange,1-fro/totalfronorm,'-o','LineWidth',2);
hold on;
plot(krange,1-RSS/TSS,'-o','LineWidth',2);
hold off;
%xlabel('k');
%ylabel('Fraction of explained variance');
hold off;
h=legend({'Frobenius','RSS'});
h.Location='southeast';
xlabel_handle.FontSize = 18;
ylabel_handle.FontSize = 18;


%% run NMF with the k value
k_final=26;
[W,H]=nnmf(V,k_final,'replicates',10);
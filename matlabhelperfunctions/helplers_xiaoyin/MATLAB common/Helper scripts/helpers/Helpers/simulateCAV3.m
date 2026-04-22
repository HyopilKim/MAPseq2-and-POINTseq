%simultate bulk injection of CAV-cre 
function [output]=simulateCAV3(inputmatrix,suppleinput,inputnormalizedmatrix,infectionthreshold,injregions,regions,plotting)


density=zeros(length(injregions),length(regions));
for i=1:length(injregions)
    infected=suppleinput(:,i)>infectionthreshold;
    density(i,:)=sum(inputnormalizedmatrix(infected,regions),1);
end
output=density./repmat(sum(density,2),1,length(regions)); %rows are injections, columns are target areas

if plotting==1
figure;
bar(output');    
xlabel('areas')
ylabel('fraction of total projection')  
end
end
    
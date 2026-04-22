function numtar_inj(barcodematrix,sourcesite,projsite,projthresh)

if exist('projthresh','var')
    barcodematrix2=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh,:);
else
    barcodematrix2=barcodematrix;
end


binary=barcodematrix2>0;

sources = barcodematrix2(:,sourcesite);          
numtar = sum(binary(:, projsite),2);
mean_tarvalue=mean(barcodematrix2(:,projsite),2);  

xvals = unique(sources);
ymean = zeros(size(xvals));
yerr = zeros(size(xvals));

for i = 1:length(xvals)
    idx = sources == xvals(i);
    y = numtar(idx);
    ymean(i) = mean(y);
    yerr(i) = std(y) / sqrt(length(y));  % standard error
end

% 에러바 그래프 그리기
figure;
errorbar(xvals, ymean, yerr, '-o', 'LineWidth', 1.5);
xlabel('Source');
ylabel('Mean number of targets');
title('Mean ± SE by source');
grid on;


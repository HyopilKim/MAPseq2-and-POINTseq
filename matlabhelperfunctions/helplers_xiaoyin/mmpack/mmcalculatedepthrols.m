function [depths,angle]=mmcalculatedepthrols(contours,lroi1,pixelsize)
% calculate distance from ROIs to both surfaces of the cortex. The rois are
% specified by a segmentation image segim1 using labels in cellid1. when
% using images of multiple resolutions, make sure that segim1 and contours
% are generated from the same file with the matching pixelsize. This
% version operates on one segim1 and one contours at a time.


ftlroi=lroi1;

[p1,d1,~]=distance2curve(contours.top,ftlroi);
[p2,d2,~]=distance2curve(contours.bottom,ftlroi);
depths=[d1 d2]*pixelsize;
u=[p1-ftlroi,zeros(size(p1,1),1)];
v=[p2-ftlroi,zeros(size(p2,1),1)];
angle=zeros(size(u,1),1);
for i=1:size(u,1)
    angle(i)=atan2d(norm(cross(u(i,:),v(i,:))),dot(u(i,:),v(i,:)));
end

save('depths.mat','depths','angle');
end
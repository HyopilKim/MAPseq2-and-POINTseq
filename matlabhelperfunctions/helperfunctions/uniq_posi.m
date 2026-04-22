function [y]=uniq_posi(x)
i=1;
y=[];
while i<length(x)
    if x(i,:)==x(i+1,:)
        y=[y;i];
        i=i+1;
    else
        i=i+1;
    end
end
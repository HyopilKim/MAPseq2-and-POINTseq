function givenproj(B,given,area)

N_area_given=sum((max(B(:,given),[],2)>0 & B(:,area)>0)>0);
N_given=sum(max(B(:,given),[],2)>0);
ratio=N_area_given/N_given;

assignin('base', 'ratio', ratio);
assignin('base', 'N_given',N_given);
assignin('base', 'N_area_given',N_area_given);

N_area_notgiven=sum((max(B(:,given),[],2)==0 & B(:,area)>0)>0);
N_notgiven=sum(max(B(:,given),[],2)==0);
ratio_not=N_area_notgiven/N_notgiven;

assignin('base', 'ratio_not', ratio_not);
assignin('base', 'N_notgiven',N_notgiven);
assignin('base', 'N_area_notgiven',N_area_notgiven);
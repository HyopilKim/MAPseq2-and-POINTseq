function [sourceratio,maxtarratio,maxs]=findinj(Bnorm,source,target,sourcetargetspikeinfold)

% input Bnorm has two columns to compare
Bnorm_source=sourcetargetspikeinfold*Bnorm(:,source);
Bnorm_target=Bnorm(:,target);
max_source=max(Bnorm_source,[],2);
min_source=min(Bnorm_source,[],2);
sourceratio=min_source./max_source;
max_target=max(Bnorm_target,[],2);
maxs=[max_source,max_target];
[mx,idx]=max(maxs,[],2);
mn=min(maxs,[],2);
ratio=mn./mx;
maxtarratio=[ratio,idx];

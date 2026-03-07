function [sourceratio,maxtarratio,maxs]=findinj(Bnorm,source,target,sourcetargetspikeinfold)

%Bnorm: normalized matrix
%source: sourcesites
%target: targetsites
%sourcetargetspikeinfold: fold difference of spikein for targets and
%injection samples for fair comparison (sourcespikein/targetspikein)
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

%% MATLAB protion of MAPseq data preprocessing.
% run preprocessing1.sh

%% You need produceBCmat2.m, normBCmat2.m

%% now, let's figure out which sequencing reads to consider for further analysis
%% first input is sample ID, second input is prefix name
% here example file names are used, so you should change the names for your
% experiments!

rankplot(1:26,'M205_HZ');


%% go back to bash

%% run produceBCmat -- this error correct and build the matrix
%% first input is sample ID, second input is prefix name, third input is name of experiment that you want make

produceBCmat2(1:12,'M205_HZ','brain1');

%%Now you have barcodematrix and spikes files ending with .mat
%clear all remove all the variables on your workspace
% and let's open the files as below or double click

clear all
load barcodematrixM205_HZ_brain1.mat
load spikesM205_HZ_brain1.mat


%% Let's normalize and filter
% B is barcode matrix with UMI

% Bseq is barcode sequences of B which written as numbers (A, T, C, G are
% repesented as specific numbers)

% Bnorm is normalized B by spikein

% sourcethresh is set as 10 so expecting at least 10 UMi in source (cell
% body)

% projthresh is set as 5 so expecting at least 5 UMI in projections
% (target)

% _filt means filtered by cellbodyThresh which set as 200 here,
% since we don't expect >200 UMI from target areas so excluding some
% barcodes showing more than 200 UMI from any target areas

% threshold values can be changed depending on experiments

% _tar means matrix including target areas only

% [B,Bnorm,B_filt,Bseq,Bnorm_filt,Bseq_filt,B_filt_tar,Bnorm_filt_tar]=normBCmat2(barcodematrix,refbarcodes,spikes,sourcethresh,projthresh,sourcesite,projsite,cellbodyThresh)
[B1,B1norm,B1_filt,B1seq,B1norm_filt,B1seq_filt,B1_filt_tar,B1norm_filt_tar]=normBCmat2(barcodematrix,refbarcodes,spikes,10,5,1,2:26,200);



%% plot some matrices!
targets=["OB","ACB","AI","CP","MTN","BLAa","PIR","VTA","TeA","ENTI","RSP","PAG"];



clustergram(X118_filt_sort,'Standardize','Row','ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')

%try a different distance metric
clustergram(Bnorm1_filt,'Standardize','Row','ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula','RowPDist','cosine','ColumnPDist','cosine')


%normalize cells by maximum projection, rather than total projection
%strength
clustergram(Bnorm1_filt./repmat(max(Bnorm1_filt,[],2),1,size(Bnorm1_filt,2)),'Standardize',0,'ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')

%plot brain 2
clustergram(Bnorm2_filt(randi(size(Bnorm2_filt,1),1000,1),:),'Standardize','Row','ColumnLabels',targets,'Cluster','column','Symmetric','false',...
    'Colormap','parula')

%% how much template switching is there? let's look at all brains together, as they were PCR'd together.
% sample 25 is the L1 control, sample 26 the h20 control


produceBCmat(1:26,'M205_HZ','all');
produceBCmatL1(1:26,'M205_HZ','all');


load barcodematrixL1M205_HZ_all.mat
load barcodematrixM205_HZ_all.mat
load spikesM205_HZ_all.mat

target_L2 = 1:24; %change target site info
target_L1 = 25; %change L1 target site info
num_spikein = zeros(1,length(spikes)); %for too many spike-in counts
for i=1:length(spikes)
    num_spikein(i) = length(spikes(i).counts2u);
end


num_target_L2 = sum(sum(barcodematrix(:,target_L2)));       % total num of L2 molecules in L2 targets
num_target_L1 = sum(sum(barcodematrixL1(:,target_L1)));         % total num of L1 molecules in L1 targets
num_spikes_L1 = sum(num_spikein(target_L1));           % total num of spike in molecules in L1 targets
num_spikes_L2 = sum(num_spikein(target_L2));           % total num of spike in molecules in L2 targets
num_templateswitching = sum(sum(barcodematrix(:,target_L1)));   % num of L2 molecules detected in L1 targets
c = num_templateswitching/(num_target_L2*(num_spikes_L1+num_target_L1));        % template switching coefficient

ratio_ts = 0.5*c*num_target_L2+c*num_spikes_L2;         % ratio of template swtiching molecules in all L2 targets


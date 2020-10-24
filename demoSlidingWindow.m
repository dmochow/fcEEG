% demo how to compute dynamic FC (sliding windows) for EEG data records

clear all; close all; clc

% download data from here
% https://www.dropbox.com/sh/ml8bsmr5v2dq1sg/AAAbdk3tqM1SOaFKbNONYQHBa?dl=0
% put it in ../data/superbowl or change the path below

pathToData='../data/superbowl/preprocessed/2012/';
filenames=dir(fullfile(pathToData,'*.mat'));
nFilenames=numel(filenames);

% parameters for FC computation
compSelect=[1 6]; % PCA, 6 components
winlenSecs=5; % sliding window length
winshiftSecs=1; % sliding window shift
fs=256; % sampling rate of the EEG
nAds=10; % # of stimuli (ads) in the demo data
gamma=0; % don't include gamma band

allFc=cell(nAds,1); % stores all subjects' FC matrices for each ad

% run sliding window fc
for f=1:nFilenames % loop across subjects
    f
    load(fullfile(pathToData,filenames(f).name));
    for a=1:nAds
        data=squeeze(Y1{a});
        data(isnan(data))=0; % no support for nans yet
        nSecs=floor(size(data,2)/fs);
        stats = fcOnly(data,compSelect,winlenSecs,winshiftSecs,fs,gamma);
        allFc{a}=cat(3,allFc{a},stats.fc);
    end
end

%%
% plot FC time series for each ad
nRows=2;
nCols=5;
figure;
for a=1:nAds
    hs(a)=subplot(nRows,nCols,a);
    muFc=nanmean(allFc{a},3); % mean across subjects
    muMuFc=nanmean(muFc,1); % mean across elements of the FC matrix
    plot(muMuFc(2:end)); % first window is artificially inflated
    axis square
    xlabel('time window');
    ylabel('mean FC');
    title(['Advertisement ' num2str(a)],'FontWeight','normal');
    
end

%% 
% show FC matrices across time for a single subject / ad
a=4; f=6; % ad 4, subject 6
fcMats=squeeze(allFc{a}(:,:,f));

hf=figure;
nRows=2; nCols=3;
times2show=2:4:24;
for t=1:numel(times2show)
    tfcMat=fcReshape(fcMats(:,times2show(t)));
    hs(t)=subplot(nRows,nCols,t);
    imagesc(tfcMat);
    title(['Time window ' num2str(times2show(t))],'FontWeight','normal');
end


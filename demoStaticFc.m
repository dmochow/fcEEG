% show how to compute a single FC matrix for EEG data records

clear all; close all; clc

% download data from here
% https://www.dropbox.com/sh/ml8bsmr5v2dq1sg/AAAbdk3tqM1SOaFKbNONYQHBa?dl=0
% put it in ../data/superbowl or change the path below

pathToData='../data/superbowl/preprocessed/2012/';
filenames=dir(fullfile(pathToData,'*.mat'));
nFilenames=numel(filenames);

% parameters for FC computation
compSelect=[1 6]; % PCA, 6 components
fs=256; % sampling rate of the EEG
nAds=10; % # of stimuli (ads) in the demo data
gamma=0; % don't include gamma band
allFc=cell(nAds,1); % stores all subjects' FC matrices for each ad

% compute FC matrix for each ad/subject
for f=1:nFilenames % loop across subjects
    f
    load(fullfile(pathToData,filenames(f).name));
    for a=1:nAds
        data=squeeze(Y1{a});
        data(isnan(data))=0; % no support for nans yet
        nSecs=floor(size(data,2)/fs);
        winlenSecs=nSecs; winshiftSecs=nSecs;
        stats = fcOnly(data,compSelect,winlenSecs,winshiftSecs,fs,gamma);
        allFc{a}=cat(2,allFc{a},squeeze(stats.fc));
    end
end

%%
% plot subject-averaged FC matrix for each ad
nRows=2;
nCols=5;
figure;
for a=1:nAds
    hs(a)=subplot(nRows,nCols,a);
    muFc=nanmean(allFc{a},2); % mean across subjects
    muFc2D=fcReshape(muFc);
    imagesc(muFc2D); % display FC matrix
    axis square
    xlabel('time window');
    ylabel('mean FC');
    title(['Advertisement ' num2str(a)],'FontWeight','normal');
    caxis([-0.5 0.5]); % show Pearson correlation from -0.5 to +0.5
end
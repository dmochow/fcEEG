function stats = fcOnly(data,compSelect,winlenSecs,winshiftSecs,fs)
%FCONLY take data and output FC matrices along time
% (c) Jacek P. Dmochowski AKA "Loquacious D", 2018-

% data: channels x time EEG array (NB: do not pass this function an EEG
%   array where the number of samples is less than the number of channels)
% compSelect: a two-element vector where the first element selects the
%   component decomposition technique (1=PCA default, 2=ICA), and the second element
%   selects the number of components (defaults to 6)
% winlenSecs: window length in seconds (defaults to 5)
% winshiftSecs: window shift in seconds (defaults to 1)
% fs: sampling rate of EEG, defaults to 64 Hz

% stats: output struct with the following fields
%   fs: sampling rate
%   fc: flattened functional connectivity matrix for each time window
%   winStarts: start times of each window

if nargin<5 || isempty(fs), fs=64; end
if nargin<4 || isempty(winshiftSecs), winshiftSecs=5; end
if nargin<3 || isempty(winlenSecs), winlenSecs=5; end
if nargin<2 || isempty(compSelect), compSelect=[1 6]; end

data = forceSpaceTime(data);
nSamples=size(data,2);

if compSelect(1)==1 % PCA
    D=compSelect(2);
    [U,S,V]=svd(data',0);
    data=U(:,1:D).';
elseif compSelect(1)==2 % ICA
    D=compSelect(2);
    [data, W, T, mu] = kICA(data,D);
else % no component analysis
    D=size(data,1);
end

%%
winlen=winlenSecs*fs;
winshift=winshiftSecs*fs;
filterOrder=fs;
bDelta = fir1(filterOrder,[1 4]/(fs/2));
bTheta = fir1(filterOrder,[4 8]/(fs/2));
bAlpha = fir1(filterOrder,[8 13]/(fs/2));
bBeta = fir1(filterOrder,[13 30]/(fs/2));

%% generate Hilbert envelopes
dataDelta = filter(bDelta,1,data,[],2);
envDelta=log(abs(hilbert(dataDelta.').'));

dataTheta = filter(bTheta,1,data,[],2);
envTheta=log(abs(hilbert(dataTheta.').'));

dataAlpha = filter(bAlpha,1,data,[],2);
envAlpha=log(abs(hilbert(dataAlpha.').'));

dataBeta = filter(bBeta,1,data,[],2);
envBeta=log(abs(hilbert(dataBeta.').'));

eData=[envDelta;...
    envTheta;...
    envAlpha;...
    envBeta];

%% windowed covariance
nWins=floor( (nSamples-winlen-1)/winshift + 1);
for w=1:nWins
    winInds=(w-1)*winshift+1:(w-1)*winshift+winlen;
    winData=eData(:,winInds);
    winR=corrcoef(winData');
    keepInds=find(triu(ones(size(winR)),1));
    fc(:,w) = winR(keepInds);
end
winStarts=(0:nWins-1)*winshift;
%

stats.fs=fs;
stats.winStarts=winStarts;
stats.fc=fc;


end


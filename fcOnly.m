function stats = fcOnly(data,compSelect,winlenSecs,winshiftSecs,fs,gamma)
%FCONLY take data and output FC matrices along time
% (c) Jacek P. Dmochowski AKA "Loquacious D" AKA "Pomeranian Boss", 2018-

if nargin<6 || isempty(gamma), gamma=0; end
if nargin<5 || isempty(fs), fs=64; end
if nargin<4 || isempty(winshiftSecs), winshiftSecs=1; end
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

if gamma && fs<110
    error('JD: sampling rate too low for gamma band');
end

%%
winlen=winlenSecs*fs;
winshift=winshiftSecs*fs;
filterOrder=fs;
bDelta = fir1(filterOrder,[1 4]/(fs/2));
bTheta = fir1(filterOrder,[4 8]/(fs/2));
bAlpha = fir1(filterOrder,[8 13]/(fs/2));
bBeta = fir1(filterOrder,[13 30]/(fs/2));
bGamma= fir1(filterOrder,[30 55]/(fs/2));

%% generate Hilbert envelopes
dataDelta = filter(bDelta,1,data,[],2);
envDelta=log(abs(hilbert(dataDelta.').'));

dataTheta = filter(bTheta,1,data,[],2);
envTheta=log(abs(hilbert(dataTheta.').'));

dataAlpha = filter(bAlpha,1,data,[],2);
envAlpha=log(abs(hilbert(dataAlpha.').'));

dataBeta = filter(bBeta,1,data,[],2);
envBeta=log(abs(hilbert(dataBeta.').'));

% 
if gamma
    dataGamma = filter(bGamma,1,data,[],2);
    envGamma=log(abs(hilbert(dataGamma.').'));
    
    eData=[envDelta;...
    envTheta;...
    envAlpha;...
    envBeta;...
    envGamma];

else
    
    eData=[envDelta;...
    envTheta;...
    envAlpha;...
    envBeta];

end



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


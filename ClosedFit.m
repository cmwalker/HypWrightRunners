function [ fits, fitErr, SNR, nLac ] = ClosedFit( base,fitParams,raw,t,freqAxis )
import HypWright.*
import HypWright.Models.*
import HypWrightRunners.*
%% Initilaize variable
if isempty(base)
    base = struct();
end
if isempty(fitParams)
    error('No initial guesses were passed in for any fit parameters');
end
Default = struct('FWHMRange', [],'endTime', 100, 'T1a', 56, 'T1b', 30,...
    'A', MultiPool(),'Kab', 0.1, 'flipAngle', 20, 'TR', 2,...
    'noiseLevel', 0, 'nAverages', 1,'lb',[],'ub',[],'verbose', false,...
    'centers',[],'fitOptions', optimset('lsqcurvefit'));
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
if isempty(base.FWHMRange) || isempty(base.centers)
    [I, ~, peakI] = FWHMRange(freqAxis, sum(abs(fftshift(fft(raw,[],1),1)),2));
    base.FWHMRange = I;
    base.centers = peakI;
end
FWHMRange = base.FWHMRange;
A = base.A;
Kab = base.Kab;
noiseLevel = base.noiseLevel;
nAverages = base.nAverages;
lb = base.lb;
ub = base.ub;
verbose = base.verbose;
centers = base.centers;
%% Fit data
fits = zeros(nAverages,length(fieldnames(fitParams)));
resnorm = zeros(nAverages,1);
for j = 1:nAverages
if(noiseLevel ~= 0)
        [noiseyData, SNR] = ApplyNoise(raw, noiseLevel,freqAxis,...
            FWHMRange,centers);
        FTData = fftshift(fft(noiseyData,[],1),1);
    else
        FTData = fftshift(fft(raw,[],1),1);
        SNR = inf;
end
peakMax = zeros(size(centers));
phases = zeros(size(centers));
signals = zeros(size(raw,2),size(centers,2));
PhaseData = zeros([size(FTData),length(centers)]);
% Phase correct and FWHM integrate each peak
%for m = 1:length(centers)
for m = 1:length(centers)
    if(centers(m) == 0)
        continue
    end
    [~,peakMax(m)] =  max(abs(FTData(centers(m),:))); % fint the maximal peak location
    phases(m) = angle(FTData(centers(m),peakMax(m))); % find the phase at the above point
    for n = 1:length(t)
        % FWHM integrate the Phased signal
        PhaseData(:,n,m) = real(exp(-1i*(phases(m)))*FTData(:,n));
    end
    [signals(:,m)] = SignalIntegration(freqAxis, squeeze(PhaseData(:,:,m)), FWHMRange(m,:));
end
nLac = sum(abs(signals(:,2)))/sum(sum(abs(signals)));
flipAngles(1,:) = base.flipAngle*pi/180;%zeros(1,size(raw,2))+base.flipAngle*pi/180;
flipAngles(2,:) = base.flipAngle*pi/180;%zeros(1,size(raw,2))+base.flipAngle*pi/180;
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
params = struct('ExchangeTerms',[0,base.Kab;0,0],'T1s',[base.T1a,base.T1b],...
    'TRList',t,'FaList',flipAngles,'fitOptions',base.fitOptions);
[fits(j,:),resultParams,allParams,resnorm(j),residual,exitflag,output,lambda,jacobian]...
    = A.fitData(params,fitParams,t,signals.',lb,ub);
end
fitErr = mean(resnorm);
if(verbose)
    % Display Model accuracy
    A.DataCompare(A,allParams,signals(1,:).',t,signals.')
    drawnow
end
end


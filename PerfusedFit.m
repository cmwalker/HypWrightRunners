function [ fits, fitErr, SNR, exitflag ] = PerfusedFit( base,fitParams,raw,t,freqAxis)
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
Default = struct('gamma', 67.262e6, 'readBandwidth', 4096, 'rfBandwidth', 5000,...
    'nPoints', 2048,'t0',0,'endTime', 100, 'T1a', 56, 'T2a', 0.02, 'T1b', 30,...
    'T2b', 0.02, 'kve', 0.02, 'vb', 0.09, 've', .91, 'ppma', -7e-6,...
    'ppmb', 7e-6,'gammaPdfA',2.8,'gammaPdfB',4.5,'scaleFactor',1,...
    'Kab', 0.1, 'flipAngle', 20, 'TR', 2,'verbose', false,...
    'FWHMRange', [], 'A', GammaBanksonModel(),'noiseLevel', 1e20,...
    'nAverages', 1, 'lb',0,'ub',100,'centers',[],'fitOptions',...
    optimset('lsqcurvefit'),'B0',3,'autoVIFNorm',true);
tmpNames = fieldnames(base);
for i = 1:numel(tmpNames) 
    if ~isfield(Default,tmpNames{i})
        fprintf('WARNING! the field %s was passesed in but does not match any of the default names. This variable wont be used!\n',tmpNames{i})
    end
end
% Fill with defaults
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
% calc FWHM range if needed
if isempty(base.FWHMRange) || isempty(base.centers)
    [I, ~, peakI] = FWHMRange(freqAxis, sum(abs(fftshift(fft(raw,[],1),1)),2));
    base.FWHMRange = I;
    base.centers = peakI;
end
FWHMRange = base.FWHMRange;
A = base.A;
bfit = @(t)padarray(gampdf(t,base.gammaPdfA,base.gammaPdfB),1,'post');
base.b = bfit;
Kab = base.Kab;
noiseLevel = base.noiseLevel;
nAverages = base.nAverages;
lb = base.lb;
ub = base.ub;
centers = base.centers;
verbose = base.verbose;
% Correct for VIF normilization if needed
if base.autoVIFNorm
    base.scaleFactor = 1.325707821912188e+03*sin(base.flipAngle)-0.258954792597810;
end
%% Fit data
fits = zeros(nAverages,length(fieldnames(fitParams)));
fitErr = zeros(nAverages,1);
for j = 1:nAverages
    if(noiseLevel ~= 0)
        [noiseyData, SNR] = ApplyNoise(raw, noiseLevel,freqAxis);
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
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
[x,~,resnorm,~,exitflag] = A.fitData(fitConstants,fitParams,t,signals.','lb',lb,'ub',ub);
resParams = fitConstants;
for i = 1:numel(tmpNames)
    resParams.(tmpNames{i}) = x(i);
end
fits(j,:) = x;
fitErr(j) = resnorm;
end
if(verbose)
    % Display Model accuracy
    figure('Position',[660 50 1040 400])
    [Y,T] = A.evaluate(resParams,linspace(0,t(end),1000),[0,0]);
    Y = Y.';
    plot(T,Y(:,1),'g',T,Y(:,2),'b',t,signals(:,1),'go',t,...
        signals(:,2),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Fit Parameters')
end
end


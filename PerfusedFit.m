function [ fits, fitErr, SNR ] = PerfusedFit( base,fitParams,raw,t )
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
    'nPoints', 2048, 'FWHMRange', 100,...
    't0',0,'endTime', 100, 'T1a', 56, 'T2a', 0.02, 'T1b', 30,...
    'T2b', 0.02, 'kve', 0.02, 'vb', 0.09, 've', .91, 'ppma', -7e-6,...
    'ppmb', 7e-6,'gammaPdfA',2.8,'gammaPdfB',4.5, 'A', BanksonModel(),...
    'Kab', 0.1, 'flipAngle', 20, 'TR', 2, 'noiseLevel', 1e20, 'nAverages', 1,...
    'lb',0,'ub',100,'verbose', false,'centers',[913,1137]);
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
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
%% Fit data
fits = zeros(nAverages,length(fieldnames(fitParams))+1);
fitErr = zeros(nAverages,1);
for j = 1:nAverages
[ noiseyData, SNR ] = ApplyNoise( raw, noiseLevel );
FTData = fftshift(fft(noiseyData,[],2),2);
% tmpSpec = abs(sum(real(fftshift(fft(raw,[],2),2))));
% [pks,centers] = findpeaks(tmpSpec);
% [~,I] = sort(pks,'descend');
% centers = centers(I(1:2));
peakMax = zeros(size(centers));
phases = zeros(size(centers));
signals = zeros(size(centers,2),size(raw,1));
% Phase correct and FWHM integrate each peak
for m = 1:length(centers)
    [~,peakMax(m)] =  max(abs(FTData(:,centers(m)))); % fint the maximal peak location
    phases(m) = angle(FTData(peakMax(m),centers(m))); % find the phase at the above point
    for n = 1:length(t)
        % FWHM integrate the Phased signal
        signals(m,n) = HypWright.fwhm(centers(m)-FWHMRange:centers(m)+FWHMRange,...
            abs(real(exp(1i*(phases(m)))*FTData(n,centers(m)-FWHMRange:centers(m)+FWHMRange))));
    end
end
signals(:,1) = 0;
tmpMax = max(max(signals));
signals = signals./tmpMax;
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
[x,~,resnorm,~] = A.fitData(fitConstants,fitParams,t,signals,'lb',lb,'ub',ub);
resParams = fitConstants;
for i = 1:numel(tmpNames)
    resParams.(tmpNames{i}) = x(i);
end
fits(j) = x(1);
fitErr(j) = resnorm*tmpMax;
end
if(verbose)
    % Display Model accuracy
    tmpFits = mean(fits(:,1));
    figure('Name',sprintf('Fit Kab: %.4f Actual Kab: %.4f',tmpFits,Kab),'NumberTitle','off',...
        'Position',[660 50 1040 400])
    [Y,T] = A.evaluate(resParams,linspace(0,t(end),1000),[0,0]);
    subplot(1,2,1)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,signals(1,:)./x(end),'go',t,...
        signals(2,:)./x(end),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Fit Parameters')
    actulParams = base;
    [Y,T] = A.evaluate(actulParams,linspace(0,t(end),1000),[0,0]);
    subplot(1,2,2)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,signals(1,:)./x(end),'go',t,...
        signals(2,:)./x(end),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Actual Parmeters')
end
end


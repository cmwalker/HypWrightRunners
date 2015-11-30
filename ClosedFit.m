function [ fits, fitErr, SNR ] = ClosedFit( base,fitParams,raw,t )
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
Default = struct('FWHMRange', 10,'endTime', 100, 'T1a', 56, 'T1b', 30,...
    'A', TwoSiteExchange(),'Kab', 0.1, 'flipAngle', 20, 'TR', 2,...
    'noiseLevel', 1e20, 'nAverages', 1,'lb',[],'ub',[],'verbose', false,...
    'centers',[763,1287]);
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
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
fits = zeros(nAverages,length(fieldnames(fitParams))+1);
resnorm = zeros(nAverages,1);
for j = 1:nAverages
[ noiseyData, SNR ] = ApplyNoise( raw, noiseLevel );
FTData = fftshift(fft(noiseyData,[],2),2);
% tmpSpec = sum(abs(fftshift(fft(raw,[],2),2)));
% [pks,centers] = findpeaks(tmpSpec);
% [~,I] = sort(pks,'descend');
% centers = centers(I(1:2));
% centers = sort(centers);
peakMax = zeros(size(centers));
phases = zeros(size(centers));
signals = zeros(size(centers,2),size(raw,1));
% Phase correct and FWHM integrate each peak
for m = 1:length(centers)
    [~,peakMax(m)] =  max(abs(FTData(:,centers(m)))); % fint the maximal peak location
    phases(m) = angle(FTData(peakMax(m),centers(m))); % find the phase at the above point
    for n = 1:length(t)
        % FWHM integrate the Phased signal
        signals(m,n) = abs(sum(real(exp(1i*(phases(m)))*FTData(n,centers(m)-FWHMRange:centers(m)+FWHMRange))));
    end
end
% tmpMax = max(max(signals));
% signals = signals./tmpMax;
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
[fits(j,:),~,resnorm(j),~] = A.fitData(fitConstants,fitParams,t,signals,lb,ub);
resParams = fitConstants;
tmp = mean(fits);
for i = 1:numel(tmpNames)
    resParams.(tmpNames{i}) = tmp(i);
end
end
fitErr = mean(resnorm);
if(verbose)
    % Display Model accuracy
    tmpFits = mean(fits(:,1));
    figure('Name',sprintf('Fit Kab: %.4f Actual Kab: %.4f',tmpFits,Kab),'NumberTitle','off',...
        'Position',[660 50 1040 400])
    [Y,T] = A.evaluate(resParams,linspace(0,t(end),1000),signals(:,1));
    subplot(1,2,1)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,signals(1,:).','go',t,...
        signals(2,:),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Fit Parameters')
    actulParams = base;
    [Y,T] = A.evaluate(actulParams,linspace(0,t(end),1000),signals(:,1).');
    subplot(1,2,2)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,signals(1,:),'go',t,...
        signals(2,:),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Actual Parmeters')
end
end


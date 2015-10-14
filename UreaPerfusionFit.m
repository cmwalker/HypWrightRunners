function [ fits, fitErr ] = UreaPerfusionFit( base,fitParams,raw,t )
import HypWright.*
import HypWright.Models.*
%% Initilaize variable
if isempty(base)
    base = struct();
end
if isempty(fitParams)
    error('No initial guesses were passed in for any fit parameters');
end
Default = struct('ureaCenter', 1025,'FWHMRange', 100,...
    't0',0,'endTime', 100, 'T1a', 56,...
     'kve', 0.06, 'vb', 0.09, 've', .91, 'ppma', 0,...
    'ppmb', 0,'gammaPdfA',2.8,'gammaPdfB',4.5, 'A', BanksonModel(),...
    'Kab', 0.0, 'flipAngle', 20, 'TR', 2, 'NoiseFact', 1e20, 'nAverages', 1,...
    'verbose', true,'lb',[],'ub',[]);
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
ureaCenter = base.ureaCenter;
FWHMRange = base.FWHMRange;
A = base.A;
bfit = @(t)padarray(gampdf(t,base.gammaPdfA,base.gammaPdfB),1,'post');
base.b = bfit;
Kab = base.Kab;
NoiseFact = base.NoiseFact;
nAverages = base.nAverages;
verbose = base.verbose;
lb = base.lb;
ub = base.ub;
%% Fit data
fits = zeros(nAverages,length(fieldnames(fitParams))+1);
fitErr = zeros(nAverages,1);
for j = 1:nAverages
noise = (rand(size(raw))-0.5)./NoiseFact;
signal = raw+noise;
FTData = fftshift(fft(signal,[],2),2);
[~,ureaPeakI] =  max(abs(FTData(:,ureaCenter)));
phaseUrea = angle(FTData(ureaPeakI,ureaCenter));
sUrea = zeros(size(t));
for i = 1:length(t)
    sUrea(i) = HypWright.fwhm(ureaCenter-FWHMRange:...
        ureaCenter+FWHMRange,abs(real(exp(1i*(phaseUrea))*...
        FTData(i,ureaCenter-FWHMRange:ureaCenter+FWHMRange))));
end
sUrea(1) = 0;
tmpMax = max(sUrea);
sUrea = sUrea./tmpMax;
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
[x,~,resnorm,~] = A.fitData(fitConstants,fitParams,t,[sUrea;zeros(size(sUrea))],...
    'lb',lb,'ub',ub);
resParams = fitConstants;
for i = 1:numel(tmpNames)
    resParams.(tmpNames{i}) = x(i);
end
fits(j,:) = x;
fitErr(j) = resnorm;
if(verbose)
    % Display Model accuracy
    figure('Name',sprintf('Fit Kab: %.4f Actual Kab: %.4f',fits(j),Kab),'NumberTitle','off',...
        'Position',[660 50 1040 400])
    [Y,T] = A.evaluate(resParams,linspace(0,t(end),1000),[0,0]);
    subplot(1,2,1)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,sUrea./x(end),'go')
    legend('Modeled Pyruvate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Fit Parameters')
    actulParams = resParams;
    actulParams.Kab = Kab;
    [Y,T] = A.evaluate(actulParams,linspace(0,t(end),1000),[0,0]);
    subplot(1,2,2)
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,sUrea./x(end),'go')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data');
    title('Actual Parmeters')
end
end
end


function [ fits, fitErr, exitflag ] = FitData( base,fitParams,signals,t )
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
Default = struct( 'T1a', 56, 'T1b', 30, 'kve', 0.02, 'vb', 0.05, 've', .95,...
    'gammaPdfA',2.8,'gammaPdfB',4.5,'t0',0,'scaleFactorGpdf',1,...
    'Kab', 0.1, 'flipAngle', 20, 'TR', 2,'verbose', false,...
    'A', GammaBanksonModel(),'nAverages', 1, 'lb',0,'ub',100,...
    'fitOptions',optimset('lsqcurvefit'));
tmpNames = fieldnames(base);
for i = 1:numel(tmpNames)
    if ~isfield(Default,tmpNames{i})
        fprintf('WARNING! the field %s was passesed in but does not match any of the default names. This variable wont be used!\n',tmpNames{i})
    end
end
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
A = base.A;
bfit = @(t)base.scaleFactorGpdf*padarray(gampdf(t-base.t0,base.gammaPdfA,base.gammaPdfB),1,'post');
base.b = bfit;
Kab = base.Kab;
lb = base.lb;
ub = base.ub;
verbose = base.verbose;
nAverages = base.nAverages;
%% Fit data
fits = zeros(nAverages,length(fieldnames(fitParams))+1);
fitErr = zeros(nAverages,1);
for j = 1:nAverages
tmpMax = max(max(signals));
signals(:,1) = [0;0];
signals = signals./tmpMax;
%% Model Results
tmpNames = fieldnames(fitParams);
fitConstants = base;
for i = 1:numel(tmpNames)
    if isfield(fitConstants,tmpNames{i})
        fitConstants = rmfield(fitConstants,tmpNames{i});
    end
end
[x,~,resnorm,~,exitflag] = A.fitData(fitConstants,fitParams,t,signals,'lb',lb,'ub',ub);
resParams = fitConstants;
for i = 1:numel(tmpNames)
    resParams.(tmpNames{i}) = x(i);
end
fits(j,:) = x;
fitErr(j) = resnorm*tmpMax;
end
if(verbose)
    % Display Model accuracy
    tmpFits = mean(fits(:,1));
    figure('Name',sprintf('Fit Kab: %.4f Actual Kab: %.4f',tmpFits,Kab),'NumberTitle','off',...
        'Position',[660 50 1040 400])
    [Y,T] = A.evaluate(resParams,linspace(0,t(end),1000),[0,0]);
    plot(T,Y(1,:),'g',T,Y(2,:),'b',t,signals(1,:)./x(end),'go',t,...
        signals(2,:)./x(end),'bo')
    legend('Modeled Pyruvate','Modeled Lactate','Simulated Pyruvate Data',...
        'Simulated Lactate Data');
    title('Fit Parameters')
end
end


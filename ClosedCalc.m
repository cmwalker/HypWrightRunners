%% Finction to run a simulation with some set of parameters and return the fit
% accuracy
function [raw,t] =  ClosedCalc(base)
import HypWright.*
import HypWright.Models.*
%% Initilaize variable
if isempty(base)
    base = struct();
end
Default = struct('gamma', 67.262e6, 'readBandwidth', 4096, 'rfBandwidth', 5000,...
    'nPoints', 2048, 'endTime', 100, 'T1a', 56, 'T2a', 0.02, 'T1b', 30,...
    'T2b', 0.02, 'ppma', -7e-6,'ppmb', 7e-6, 'Kab', 0.1, 'flipAngle', 20,...
    'TR', 2,'A', TwoSiteExchange(), 'verbose', false);
tmpNames = fieldnames(Default);
for i = 1:numel(tmpNames)
    if ~isfield(base,tmpNames{i})
        base.(tmpNames{i}) = Default.(tmpNames{i});
    end
end
base.flipAngle = base.flipAngle*pi/180;
gamma = base.gamma;
readBandwidth = base.readBandwidth;
rfBandwidth = base.rfBandwidth;
nPoints = base.nPoints;
endTime = base.endTime;
T1a = base.T1a;
T2a = base.T2a;
T1b = base.T1b;
ppma = base.ppma;
ppmb = base.ppmb;
Kab = base.Kab;
flipAngle = base.flipAngle;
TR = base.TR;
verbose = base.verbose;
%TODO add input validation;
%% Init World
world = HypWright.World.getWorld;
world.initWorld()
tmp = world.getB0;
B0 = tmp(3);
Spin = TwoSiteExchangeGroup([0;0;1;0;0;0],[0;0;0;0;0;0],...
                T1a,T2a,ppma,T1b,T2a,ppmb,gamma,1,Kab,[]);
V = Voxel([0;0;0],Spin);
world.addVoxel(V);
%% Build Pulse Sequence
PS = PulseSequence;
t = 0.01:TR:endTime;
ADC = zeros(length(t),nPoints);
for i = 1:length(t)
    Pulse = SincPulse(t(i),rfBandwidth,flipAngle/(gamma),gamma*B0,[],...
        sprintf('Excitation%d',1));
    PS.addPulse(Pulse)
    ADC(i,:) = Pulse.endTime:1/readBandwidth:Pulse.endTime+(nPoints-1)/readBandwidth;
end
world.setPulseSequence(PS)
%% Calculate
world.calculate(t(end)+10);
FID = zeros(length(t),nPoints);
for i = 1:length(t)
    [FID(i,:), freqAx] = world.evaluate(ADC(i,:),-gamma*B0);
end
raw = FID;
t = linspace(ADC(1,floor(end/2)),ADC(end,floor(end/2)),size(ADC,1));
if (verbose)
    figure
    surf(freqAx,t,abs(fftshift(fft(FID,[],2),2)));
    drawnow
end
end
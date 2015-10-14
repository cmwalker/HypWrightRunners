%% Finction to run a simulation with some set of parameters and return the fit
% accuracy
function [raw,t] =  PerfusedCalc(base)
import HypWright.*
import HypWright.Models.*
%% Initilaize variable
if isempty(base)
    base = struct();
end
Default = struct('gamma', 67.262e6, 'readBandwidth', 4096, 'rfBandwidth', 5000,...
    'nPoints', 2048, 'pyrCenter', 913, 'lacCenter', 1137, 'FWHMRange', 100,...
    't0',0,'endTime', 100, 'T1a', 56, 'T2a', 0.02, 'T1b', 30,...
    'T2b', 0.02, 'kve', 0.02, 'vb', 0.09, 've', .91, 'ppma', -7e-6,...
    'ppmb', 7e-6,'gammaPdfA',2.8,'gammaPdfB',4.5, 'A', BanksonModel(),...
    'Kab', 0.1, 'flipAngle', 20, 'TR', 2, 'NoiseFact', 1e20, 'nAverages', 1,...
    'verbose', false);
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
t0 = base.t0;
T1a = base.T1a;
T2a = base.T2a;
T1b = base.T1b;
kve = base.kve;
vb = base.vb;
ve = base.ve;
ppma = base.ppma;
ppmb = base.ppmb;
A = base.A;
bfit = @(t)padarray(gampdf(t,base.gammaPdfA,base.gammaPdfB),1,'post');
base.b = bfit;
b = @(t)padarray(padarray(gampdf(t-t0,base.gammaPdfA,base.gammaPdfB),2),1,'post');
Kab = base.Kab;
flipAngle = base.flipAngle;
TR = base.TR;
verbose = base.verbose;
Mz = [];
%TODO add input validation;
%% Init World
world = HypWright.World.getWorld;
world.initWorld()
tmp = world.getB0;
B0 = tmp(3);
Spin = TwoSitePerfusionExchangeGroup([0;0;0;0;0;0],[0;0;0;0;0;0],...
    T1a,T2a,ppma,T1b,T2a,ppmb,gamma,ve,Kab,[],kve/ve,b);
Spin2 = BanksonSpinGrp([0;0;0],[0;0;0],T1a,T2a,gamma,ppma,vb,[2.8,4.5],t0);
V = Voxel([0;0;0],Spin);
V.addSpin(Spin2);
if (verbose)
V.debug = true;
end
world.addVoxel(V);
%% Build Pulse Sequence
PS = PulseSequence;
t = 0:TR:endTime;
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
if (verbose)
V.debug = true;
tMz = 0:0.1:endTime;
world.evaluate(tMz);
load('tmp')
evSpace = [];
vSpace = [];
for i = 1:size(tmp,1)
evSpace = [evSpace,tmp{i,1}];
vSpace = [vSpace,tmp{i,2}];
end
M(1:3,:) = ve*evSpace(1:3,:)+vb*vSpace;
M(4:6,:) = ve*evSpace(4:6,:);
MzPyr = M(3,:);
MzLac = M(6,:);
V.debug = false;
end
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
    figure('Name',sprintf('Kab: %.4f Flip Angle %2f Repetition Time %.4f',...
    Kab,flipAngle,TR),'NumberTitle','off','Position',[660 50 1040 400])
    plot(tMz,MzPyr,'go',tMz,MzLac,'bo')
    legend('Simulated Pyruvate Mz','Simulated Lactate Mz');
        % Display model and Mz
    figure('Name',sprintf('Kab: %.4f Flip Angle %2f Repetition Time %.4f',...
    Kab,flipAngle,TR),'NumberTitle','off','Position',[660 50 1040 400])
    [Y,T] = A.evaluate(base,linspace(0,t(end),1000),[0,0]);
    normFact = max([max(MzPyr) max(MzLac)])/max(Y(:));
    plot(T,Y(1,:),'g',T,Y(2,:),'b',tMz,MzPyr./normFact,'go',tMz,...
        MzLac./normFact,'bo')
    legend('Modeled Pyruvate','Modeled Lactate',...
        'Simulated Pyruvate Mz','Simulated Lactate Mz');
    title('Fit Parameters')
end
end
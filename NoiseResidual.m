function [ noiseResidual ] = NoiseResidual( noiseLevel, nSamples )
%NOISERESIDUAL Summary of this function goes here
%   Detailed explanation goes here
import HypWrightRunners.*
import HypWright.*
world = HypWright.World.getWorld;
world.initWorld()
tmp = world.getB0;
B0 = tmp(3);
Spin = IsolatedSpinGrp([0;0;1],[0;0;1],500,1,67.262e6,0,1);
V = Voxel([0;0;0],Spin);
world.addVoxel(V);
%% Build Pulse Sequence
PS = PulseSequence;
Pulse = SincPulse(0.1,5000,90/(67.262e6),67.262e6*B0,[],...
    sprintf('Excitation%d',1));
PS.addPulse(Pulse)
ADC = Pulse.endTime:1/4096:Pulse.endTime+(2048-1)/4096;
world.setPulseSequence(PS)
%% Calculate
world.calculate(10);
[FID, freqAxis] = world.evaluate(ADC,-67.262e6*B0);
raw = FID;
tmpFTData = fftshift(fft(raw,[],1),1);
[I, ~, peakI] = FWHMRange(freqAxis, abs(tmpFTData));
phases = angle(tmpFTData(peakI));
signals = zeros(1,nSamples);
noiseResidual = zeros(1,10);
for j = 1:10
    for i = 1:nSamples
        [noiseyData, ~] = ApplyNoise(raw, noiseLevel, freqAxis);
        tmpFTData = fftshift(fft(noiseyData,[],1),1);
        signals(i) = SignalIntegration(freqAxis,...
            real(exp(-1i*(phases))*tmpFTData), I);
    end
    fun = @(x,t)zeros(size(t))+x;
    opts = optimset('Display','off');
    [~,noiseResidual(j),~] = lsqcurvefit(fun,2,1:nSamples,signals,[],[],opts);
end
noiseResidual = mean(noiseResidual);
end


function [ noiseyData, SNR ] = ApplyNoise( raw, noiseLevel )
%APPLYNOISE Applys some level of noise to a clean dataset
%   [ noiseyData, SNR ] = ApplyNoise( raw, noiseLevel ) - applies noise scaled by
%   some factor noiseLevel to some set of raw data. SNR is the theroretical SNR 
%   for that noise factor to determin how noiseLevel will affect SNR see
%   DynamicSNR
import HypWrightRunners.*
noise = (rand(size(raw))-0.5+1i*rand(size(raw))-0.5)/noiseLevel;
noiseyData = raw+noise;
SNR = DynamicSNR(raw,noiseLevel);
end


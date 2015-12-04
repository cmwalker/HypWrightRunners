function [ SNR,meanSNR ] = DynamicSNR(raw,noise,freqAxis,varargin)
%DYNAMICSNR gives the SNR and mean SNR of a dynamic hyperpolarized study
%   [ SNR,meanSNR ] = DynamicSNR(raw,t,noiseLevel) - takes some raw noisless
%   FIDs, the time vector associated with the FIDs, and some noise factor and
%   returns the total SNR and the average SNR per fid. raw should be a m x n
%   matriz where m is the number of fids and n is the number of points per fid,
%   t should be some time vector associated witht the FID and should be m x 1
%   (actually as long as t is a vector of the size m x 1 this function will run
%   fine). Noise level is some consatant factor that applies to the noise added
%   to the raw data.
% import libraries
import HypWright.*
import HypWrightRunners.*
% Find Peak Centers
% generate noise
FTData = fftshift(fft(raw,[],1),1); % Fourier transfor data
if nargin <= 4
if isvector(FTData)
    [I, ~, peakI] = FWHMRange(freqAxis, abs(FTData));
else
    [I, ~, peakI] = FWHMRange(freqAxis, sum(abs(FTData),2));
end
else
    I = varargin{1};
    peakI = varargin{2};
end
peakMax = zeros(size(peakI));
phases = zeros(size(peakI));
signals = zeros(size(raw,2),size(peakI,2));
PhaseData = zeros([size(FTData),length(peakI)]);
for m = 1:length(peakI)
    if(peakI(m) == 0)
        continue
    end
    [~,peakMax(m)] =  max(abs(FTData(peakI(m),:))); % fint the maximal peak location
    phases(m) = angle(FTData(peakI(m),peakMax(m))); % find the phase at the above point
    for n = 1:size(FTData,2)
        % FWHM integrate the Phased signal
        PhaseData(:,n,m) = real(exp(-1i*(phases(m)))*FTData(:,n));
    end
    [signals(:,m)] = SignalIntegration(freqAxis, squeeze(PhaseData(:,:,m)), I(m,:));
end
totalSignal = sum(sum(signals)); % Find total raw signal
noiseFT = fftshift(fft(noise,[],1),1); % find FFT of raw noise
% estimate the avearg std of the noise
noiseEstimate = mean(std(noiseFT,[],1)); 
SNR = totalSignal/noiseEstimate; % calculate total SNR
meanSNR = SNR/size(raw,2); % calculate average SNR per timepoint
end


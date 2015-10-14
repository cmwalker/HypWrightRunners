function [ SNR,meanSNR ] = DynamicSNR(raw,noiseLevel)
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
% Find Peak Centers
if size(raw,1) == 1
    [~,centers] = findpeaks(abs(fftshift(fft(raw,[],2),2)));
else
    [~,centers] = findpeaks(abs(sum(fftshift(fft(raw,[],2),2))));
end
% generate noise
noise = (rand(size(raw))+1i*rand(size(raw)))/noiseLevel;
FTData = fftshift(fft(raw,[],2),2); % Fourier transfor data
t = 1:size(raw,1);
% Initialize variables 
peakMax = zeros(size(centers));
phases = zeros(size(centers));
signals = zeros(size(centers,2),size(raw,1));
% Phase correct and FWHM integrate each peak
for i = 1:length(centers)
    [~,peakMax(i)] =  max(abs(FTData(:,centers(i)))); % fint the maximal peak location
    phases(i) = angle(FTData(peakMax(i),centers(i))); % find the phase at the above point
    for j = 1:length(t)
        % FWHM integrate the Phased signal
        signals(i,j) = HypWright.fwhm(centers(i)-100:centers(i)+100,...
            real(exp(1i*(phases(i)))*FTData(j,centers(i)-100:centers(i)+100)));
    end
end
totalSignal = sum(sum(signals)); % Find total raw signal
noiseFT = fftshift(fft(noise,[],2),2); % find FFT of raw noise
% estimate the avearg std of the noise
noiseEstimate = mean(std(noiseFT(:,round(1:end/3)),[],2)); 
SNR = totalSignal/noiseEstimate; % calculate total SNR
meanSNR = SNR/length(t); % calculate average SNR per timepoint
end


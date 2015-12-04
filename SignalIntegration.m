function [ signals ] = SignalIntegration(xAxis, rawSignal, rangeIndices)
%SIGNALINTEGRATION Summary of this function goes here
%   Detailed explanation goes here
signals = zeros(size(rawSignal,2),size(rangeIndices,1));
for i = 1:(size(rangeIndices,1))
    for j = 1:size(rawSignal,2)
        tmpSpec = interp1(1:length(xAxis),rawSignal(:,j),...
            [rangeIndices(i,1),ceil(rangeIndices(i,1)):floor(rangeIndices(i,2)),rangeIndices(i,2)]);
        tmpPpm = interp1(1:length(xAxis),xAxis,...
            [rangeIndices(i,1),ceil(rangeIndices(i,1)):floor(rangeIndices(i,2)),rangeIndices(i,2)]);
        signals(j,i) = trapz(tmpPpm,tmpSpec);
    end
end
end


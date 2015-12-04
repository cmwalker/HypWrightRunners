function [I, C, peakI] = FWHMRange(xAxis, yAxis, varargin)
%FWHMRANGE Summary of this function goes here
%   Detailed explanation goes here

% parse input
p = inputParser;
addOptional(p,'verbose',false);
parse(p,varargin{:});
verbose = p.Results.verbose;

% Find peaks
peaks = findpeaks(yAxis);
halfHieghts = peaks./2;
peakI = zeros(length(peaks),1);
for i = 1:length(peaks)
    peakI(i) = find(yAxis==peaks(i));
end
% Split area between peaks
breaks = zeros(length(peaks)-1,1);
for i = 1:length(breaks)
    breaks(i) = floor(mean(peakI(i:i+1)));
end
breaks = [1,breaks,length(yAxis)]; % pad break pints
if verbose
    figure('Name','FWHM Ranges','Position',[440,140,760,660])
    ax1 = subplot(2,1,1);
    plot(yAxis)
    xlabel('Indices'),ylabel('Signal'),title('Indexed')
    hold(ax1,'on')
    ax2 = subplot(2,1,2);
    plot(xAxis,yAxis)
    xlabel('x Values'),ylabel('Signal'),title('User supplied X-Axis')
    hold(ax2,'on')
end
% loop trough the peak segments
I = zeros(length(peaks),2);
C = zeros(length(peaks),2);
for i = 1:(length(peaks))
    right = find(yAxis(breaks(i):breaks(i+1))>=halfHieghts(i),1,'first');
    right = right+breaks(i)-1;
    righti = interp1([yAxis(right-1),yAxis(right)],[right-1,right],halfHieghts(i));
    left = find(yAxis(breaks(i):breaks(i+1))>=halfHieghts(i),1,'last');
    left = left+breaks(i)-1;
    lefti = interp1([yAxis(left),yAxis(left+1)],[left,left+1],halfHieghts(i));
    I(i,:) = [righti,lefti];
    C(i,:) = interp1(1:length(xAxis),xAxis,[righti,lefti]);
    if verbose
        plot(ax1,righti,halfHieghts(i),'ro')
        plot(ax1,lefti,halfHieghts(i),'ko')
        plot(ax2,C(i,1),halfHieghts(i),'ro')
        plot(ax2,C(i,2),halfHieghts(i),'ko')
    end;
end
if verbose
    legend(ax2,'Signal','Left FWHM Ranges','Right FWHM Ranges');
    legend(ax1,'Signal','Left FWHM Ranges','Right FWHM Ranges');
    hold(ax1,'off')
    hold(ax2,'off')
end
end


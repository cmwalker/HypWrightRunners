function [fitKpl, fitErr, raw, t] = PerfusedAccuracy(base,fitParams )
%PERFUSEDACCURACY2 Summary of this function goes here
%   Detailed explanation goes here
import HypWrightRunners.*
[raw,t] =  PerfusedCalc(base);
[fitKpl, fitErr] = PerfusedFit( base,fitParams,raw,t);
end


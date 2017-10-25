% Estimate location of <x> using Huber's M-estimator with parameter 
% <huberParam>. Return NaN if any element in <x> is a NaN.
%
% x: Input data as a numeric vector
% huberParam: Parameter of Huber's M-estimator
%_______________________________________________________________________
% Copyright (C) 2013 Camille Maumet
function y = robust_exclude_nan(x, huberParam)
    if any(isnan(x)) 
        y = NaN;
    else
        y = robustfit(ones([numel(x) 1]), x, 'huber', huberParam, 'off');
    end
end
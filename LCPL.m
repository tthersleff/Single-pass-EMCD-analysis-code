function [params, BG, BGsub] = LCPL(energy, spectrum, window, varargin)
%LCPL Linear Combination of Power Laws
%   LCPL returns the fitting parameters for a linear fit of power laws.
%   This is useful for modelling of the pre-edge EELS spectrum.  Inspired
%   by the paper of Cuevera et al. detailing the Cornell Spectrum Imager
%   software and is adopted to Matlab here
%
%   (c) 2019 Thomas Thersleff, Stockholm University

r = [-3 -2];
method = 'fast'; %'fast', 'constrained', 'robust'

% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'r'
                r = varargin{idx + 1};
            case 'method'
                method = varargin{idx + 1};
        end
    end
end

preE = energy(makeWin(energy, window));
preS = spectrum(makeWin(energy, window));

X(:, 1) = preE .^ r(1);
X(:, 2) = preE .^ r(2);
if strcmpi(method, 'fast')
    params.c = X\preS;  
elseif strcmpi(method, 'robust')
    error('Robust fitting method currently not implemented');
    % params.c = robustfit(X, preS, 'bisquare', 'const', 'off');
    % params.c = robustfit(X, preS);
elseif strcmpi(method, 'constrained')
    % This is useful if you want to constrain the fitting parameters to be
    % non-negative.  It increases the processing time considerably
    % (approximately 5 times longer, by my testing)
    opts = optimoptions('lsqlin', 'Display', 'off');
    params.c = lsqlin(X, preS, [], [], [], [], zeros(1, 2), [inf inf], [], opts);
else
    error('Unknown regression method');
end

params.r = r;

if nargout >= 2
    X2(:, 1) = energy .^ r(1);
    X2(:, 2) = energy .^ r(2);
    BG = params.c(1:2)' * X2';
    BG = BG(:);
end

if nargout == 3
    BGsub = spectrum - BG;
end

end


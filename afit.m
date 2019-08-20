function [coeffs, BG, BGsub] = afit(energy, spectrum, window, varargin)
%AFIT Simple power-law fit for EELS data
%   afit returns the parameters of a simple power-law fit to an EELS
%   spectrum.  Based off of Egerton's Afit routine.
%
%   (c) 2019 Thomas Thersleff, Stockholm University

method = 'fast'; %'fast', 'robust'

% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'method'
                method = varargin{idx + 1};
        end
    end
end

% Crop the spectra
preE = energy(makeWin(energy, window));
preS = spectrum(makeWin(energy, window));

% % Take the log
% Elog = log(preE);
% Slog = log(preS);

% Fit a line
if strcmpi(method, 'fast')
    coeffs = [ones(size(preE)) log(preE)] \ log(preS);
elseif strcmpi(method, 'robust')
    coeffs = robustfit(log(preE), log(preS));
end

% Adjust the parameters
coeffs(1) = exp(coeffs(1));

if nargout >= 2
    BG = coeffs(1) .* energy .^ coeffs(2);
end

if nargout == 3
    BGsub = spectrum - BG;
end

end


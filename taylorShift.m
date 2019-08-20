function [Z, shifts] = taylorShift(Y, varargin)
%TAYLORSHIFT Calculates the shift of a signal using Taylor expansion
%regression
%   taylorShift calculates the shift of a signal by computing the Taylor
%   coefficients of a regression between the signal mean and its first
%   derivative.  It is based off of the work by Witjes et. al:

%   H. Witjes, M. van den Brink, W. J. Melssen, and L. M. C. Buydens, 
%   Chemometrics and Intelligent Laboratory Systems 52, 105 (2000).

%   This version is quite primative.  The input Y must be a 2D matrix.  If
%   the signal of interest has a low SNR, then it results in large
%   outliers, which cause problems.  So a mask can be included in the
%   input arguments to avoid this.

%   (c) 2019 Thomas Thersleff, Stockholm University


%% Preliminary

% Set default values
nIter = 10;
mask = true(1, size(Y, 2));
order = 1;
robust = false;
snrThresh = [];

% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'mask'
                mask = varargin{idx + 1};
            case 'niter'
                nIter = varargin{idx + 1};
            case 'order'
                order = varargin{idx + 1};
            case 'robust'
                robust = varargin{idx + 1};
            case 'snr'
                snrThresh = varargin{idx + 1};
        end
    end
end

if order > 1
    fprintf('\nCorrection for higher order than peak shift.  Number of iterations constrained to 1\n');
    nIter = 1;
end

%% Perform the iteration

Z = Y;
% if ~exist('mask', 'var')
%     mask = true(1, size(Y, 2));
% end

if ~isempty(snrThresh)
    mu = smooth(mean(Z(:, mask)', 'omitnan'));
    projMu = mu\Z; % Project the mean onto the raw data
    scaledMu = mu * projMu; % Scale the mean to fit the raw data
    ns = Z - scaledMu; % Get the noise matrix
    nPix = size(Z, 2); % Number of "pixels"
    for ind = nPix:-1:1
        SP(ind) = norm(Z(:, ind)); % Get Signal Power
        NP(ind) = norm(ns(:, ind)); % Get Noise Power
    end
    SNRmap = 10*log10(SP ./ NP); % Get the SNR map

%     shifts(SNRmap < snrThresh, idx) = 0; % Remove values lower than SNRthresh
    mask = SNRmap > snrThresh;
end



for idx = 1:nIter
    
    fprintf(1, '\nIteration number %02d', idx);
    X(:, 1) = mean(Z(:, mask)');
    mu = smooth(mean(Z(:, mask)', 'omitnan'));
    D = mu;
%     X(:, 1) = smooth(mean(Z', 'omitnan'));
    for ind = 1:order
        D = gradient(D);
    end
    
    X(:, 1) = mu;
    X(:, 2) = D;
%     X(:, 3) = gradient(X(:, 2));

%     B = inv(X' * X) * X' * Z;
    B = X\Z;
    
    shifts(:, idx) = (B(2, :) ./ B(1, :))';
%     broadening(:, idx) = (B(3, :) ./ B(1, :))';
    
%     bad = isoutlier(shifts, 'movmean', [5 5]);
%     bad = abs(shifts) > maxShift;
%     mm = movmedian(shifts, [5 5]);
%     shifts(bad) = mm(bad);
    
%     shifts(isnan(shifts)) = 0;
    
%     G = X(:, 3) * broadening(:, idx)';
%     Z = Z - G;
    
%     if robust
%         shifts(~mask, idx) = NaN;
%         shifts(:, idx) = filloutliers(shifts(:, idx), 'linear', 'median');
%         shifts(~mask, idx) = 0;
%     else
%         shifts(~mask, idx) = 0;
%     end

    shifts(~mask, idx) = 0;
    
    isomap.shift = sum(shifts, 2);
    
    if order == 1
        Z = applyISOmap(Y, isomap, 'crop', true);
    end
    
%     Z(isnan(Z)) = 0;
    
    clear X

end


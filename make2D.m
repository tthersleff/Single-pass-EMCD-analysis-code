function [ ESI2D, reshapeStruct ] = make2D( ESI, signalDims , varargin)
%GET2DESI Makes a 2D ESI, useful for many Matlab functions
%   get2DESI makes any ESI into a 2D matrix, which is useful for many
%   Matlab functions and for loops.  Importantly, it does not throw an
%   error if a 1D vector is passed, so it can be easily integrated into all
%   functions
%   (c) 2017 Thomas Thersleff, Stockholm University
%
%   ESI --> The ESI to be turned into a 2D matrix.  Can be any size (1D,
%   2D, 3D, 4D, etc).
%   signalDims (optional) --> The dimensions of the signal axes.  Default:
%   1
%   transposeFlag (optional) --> transpose the 2D matrix to have dimensions
%   [navigationAxis signalAxis].  Default: 0 (false)
%   ESI2D --> The 2D ESI.  Dimensions [signalAxis navigationAxis]
%   reshapeStruct --> a struct that can be passed to the makeTensor
%   function to recover the original shape and dimensions of the ESI.

%% Check inputs

% Optional input signalDims
if ~exist('signalDims', 'var')
    signalDims = 1;
end

% Set defaults for argument pair inputs
transposeFlag = 0;
navDims = []; % Empty to trigger the if statement if not explicitly passed

% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'transpose'
                transposeFlag = varargin{idx + 1};
            case 'navfirst'
                transposeFlag = varargin{idx + 1};
            case 'navdims'
                navDims = varargin{idx + 1};
        end
    end
end

%%  Prepare for reshape

dims = size(ESI);
allDims = 1:ndims(ESI);

% Distinguish between signal and navigation axes:
isSignal = false(1, ndims(ESI));
isSignal(signalDims) = true;

% If no navigation axes were explicitly passed, they will be set to all
% non-signal axes
if isempty(navDims)
    navDims = allDims(~isSignal);
end

% Get the reshape size
nSignal = prod(dims(isSignal));
nNav = prod(dims(~isSignal));
reshapeOrd = [nSignal nNav];

% Set the permutation order
if transposeFlag
    permOrd = [navDims signalDims];
    reshapeOrd = fliplr(reshapeOrd);
else
    permOrd = [signalDims navDims];
end

%% Permute and reshape

ESI2D = permute(ESI, permOrd);
ESI2D = reshape(ESI2D, reshapeOrd);

if nargout == 2
    reshapeStruct.dims = dims;
    reshapeStruct.permOrd = permOrd;
    reshapeStruct.reshapeOrd = reshapeOrd;
    reshapeStruct.transposeFlag = transposeFlag;
    reshapeStruct.allDims = allDims;
    reshapeStruct.isSignal = isSignal;
    reshapeStruct.nPix = nNav;
    reshapeStruct.nE = nSignal;
    reshapeStruct.navDims = dims(navDims);
end


end


function spikeMap = getXRaySpikeModel(ESI, varargin)
%GETXRAYSPIKEMODEL Identify and remove X-Ray spikes in EELS datacubes
%   getXRaySpikeModel is a function intended to identify and remove x-ray
%   spikes in EELS spectrum images.  These spikes are identified as
%   outliers based on computing a moving median of the entire dataset.
%   This version iterates along the spectral dimension, meaning that the
%   x-rays are identified in the spatial dimensions.  Note that the image
%   processing toolbox is required to perform dialation.  And tensor lab is
%   used for the sparse representation.  If these are not installed, I have
%   some minimal error checking to ensure that the function doesn't break.

%   Tensorlab can be downloaded from:
%   https://www.tensorlab.net/

%   INPUTS
%   ESI:    3D EELS spectrum image, dimensions [E x y] (energy axis first)

%   OUTPUTS
%   spikeMap:   Sparse representation of the x-ray spikes

%   (c) 2018 Thomas Thersleff, Stockholm University

%% Set default parameters

% Default parameters
nSigma = 10;
medWin = 5;
showPlots = false;
method = 'replace';
dilationLength = 1;

% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'sigma'
                nSigma = varargin{idx + 1};
            case 'medwin'
                medWin = varargin{idx + 1};
            case 'plot'
                showPlots = varargin{idx + 1};
            case 'method'
                method = varargin{idx + 1};
            case 'dilate'
                dilationLength = varargin{idx + 1};
        end
    end
end

tensorlab;

%% Get moving median (spike replacements)

med = movmedian(ESI, medWin, 1, 'Endpoints', 'shrink'); %Median ESI

%% Prepare for spike identification

% Subtract out median from raw
sub = abs(ESI-med);

% Preallocation for speed
nChannels = size(ESI, 1);
stddev = zeros(nChannels, 1); 

% Get standard deviation in spatial mode
for ind = 1:nChannels
    stddev(ind, 1) = std2(ESI(ind, :, :));
end

%% Find spikes
%This will give an error on Matlab versions less than 2016!
isSpike = sub >= nSigma.*stddev; %Find the spikes

% % Get linear indices and subscripts of the spikes
% idxSpike = find(isSpike);
% [iE, iRow, iCol] = ind2sub(size(isSpike), idxSpike);
% 
% %Dilate the NaN mask surrounding the spike using a morphological operator
% se = strel([1;1;1]); %Structure element, column wise
% isSpike2D = imdilate(make2D(isSpike), se); 

%% Dilate the spike matrix

if dilationLength > 1
    try
        se = strel('cuboid', [dilationLength 1 1]);
        isSpike = imdilate(isSpike, se);
    catch
        fprintf(1, '\nThe image processing toolbox in Matlab is required for the dialate function to work.  If this is not installed, the logical spike mask will not be dialated');
    end
end
    
%% Replace the spikes

switch lower(method)
    case 'replace'
        ESIout = ESI;
        ESIout(isSpike) = med(isSpike);
    case 'interpolate'
        % Replace spikes with NaN values
        ESInan = ESI;
        ESInan(isSpike) = NaN;
        
        % Paint over the NaN values
        [D, rs] = make2D(ESInan);
        ESIout = inpaint_nans(D);
        ESIout = reshapeESI(ESIout, rs);
end


%% Get the spike map

try
    spikeMap = fmt(ESI - ESIout); % Sparse matrix, tensorlab implementation
catch
    % If tensorlab is not installed, use Matlab's built-in sparse matrix
    % representation
    spikeMap = sparse(ESI - ESIout);
end

%% Plot if desired
try
    if showPlots
        try
            viewESI(1:size(ESI, 1), 1, ful(spikeMap), ESIout);
        catch
            % If tensorlab is not installed, use Matlab's built-in converter
            % from sparse to full matrices
            viewESI(1:size(ESI, 1), 1, full(spikeMap), ESIout);
        end
    end
catch
    % The GUI Layout toolbox can be installed from:
    % http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
    fprintf(1, 'The viewESI function requires the GUI Layout toolbox');
    showPlots = 0;
end

end


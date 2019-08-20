function isomap = getISOmap(ESI, varargin)
%GETISOMAP Produces an isochromaticity map for EELS spectrum images
%   getISOmap produces a structure (isomap) that contains information about
%   the energy drift offset for a hyperspectral EELS datacube.  The energy
%   drift offset is determined by cross correlating a reference spectrum
%   (ref --> optional input) to all of the EELS spectra in the datacube.
%   It thus closely mimicks the behavior of Digital Micrograph.  The
%   generated map is applied using the function applyISOmap.m

%   Note that, if you want to use the high pass filter, then you need both
%   the DSP and signal processing toolboxes installed.  I apologize if you
%   don't have these.  If this is a problem, either write some code to do
%   it (and let me know if you do!) or simply set highPassWidth = 0.

%   (c) 2018 Thomas Thersleff, Stockholm University


%% Import the data

% Set default values
signalDims = 1;
smoothData = true;
ref = 1;
transpose = false;
diffOrder = 0;
interp = 1;
highPassWidth = 0;
signalMask = true(size(ESI, 1), 1);
showProgress = false;
useParallel = false;


% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'signaldims'
                signalDims = varargin{idx + 1};
            case 'smooth'
                smoothData = varargin{idx + 1};
            case 'ref'
                ref = varargin{idx + 1};
            case 'transpose'
                transpose = varargin{idx + 1};
            case 'difforder'
                diffOrder = varargin{idx + 1};
            case 'interp'
                interp = varargin{idx + 1};
            case 'highpasswidth'
                highPassWidth = varargin{idx + 1};
            case 'signalmask'
                signalMask = varargin{idx + 1};
            case 'progressbar'
                showProgress = varargin{idx + 1};
            case 'parallel'
                useParallel = varargin{idx + 1};
        end
    end
end

% Get the ESI
[Do, rs] = make2D(ESI(signalMask, :, :), signalDims, 'transpose', transpose);
D = Do;

% Get the reference. If it's a linear index, do nothing
if length(ref) == 1
    refidx = ref;
    R = D(:, refidx);
elseif length(ref) == 2
    refidx = sub2ind(rs.navDims, ref(1), ref(2));
    R = D(:, refidx);
else
    R = ref;
end



%% Prepare the data

if highPassWidth >= 2
    [b, a] = butter(1, 1/highPassWidth, 'high');
    R = filtfilt(b, a, R);
    for idx = 1:rs.nPix
        D(:, idx) = filtfilt(b, a, D(:, idx));
    end
end

if smoothData
    R = smooth(R);
    for idx = 1:rs.nPix
        D(:, idx) = smooth(D(:, idx));
    end
end

if diffOrder ~= 0
    R = diff(R);
    D = diff(D, diffOrder, 1);
end

%% Perform cross correlation

% shiftRaw = zeros(rs.nPix, 1);
% coeff = zeros(rs.nPix, 1);

% for idx = 1:rs.nPix
%     xc = xcorr(R, D(:, idx), 'coeff');
%     if interp >= 2
%         xc = resampleESI(xc, [length(xc)*interp, 1], 'spline');
%     end
%     
%     [coeff(idx), M] = max(abs(xc));
%     shiftRaw(idx) = (M - length(R)*interp)/interp;
% end

    
dftresult = zeros(rs.nPix, 4);

% if showProgress; fprintf(1, 'Progress:   0%%'); end
tic;
if useParallel
    fprintf(1, '\nGetting isomap using parallel computing toolbox...');
    if showProgress
        PP = parallel.pool.DataQueue;
        h = waitbar(0, 'Progress');
        afterEach(PP, @nUpdateWaitbar);
        pp = 1;
    end
    
    parfor idx = 1:rs.nPix
        dftresult(idx, :) = dftregistration(fft(R), fft(D(:, idx)), interp);
%         if showProgress; send(PP, idx); end
        if showProgress; fprintf(1, '\n%06d', idx); end
    end
    
    if showProgress; delete(h); end
    
else
    fprintf(1, '\nGetting isomap without parallel computing toolbox...');
    if showProgress; fprintf(1, 'Progress:   0%%'); end
    for idx = 1:rs.nPix
        dftresult(idx, :) = dftregistration(fft(R), fft(D(:, idx)), interp);
        if showProgress
            fprintf(1, '\b\b\b\b\b%4d%%', round((idx/rs.nPix) * 100));
        end
    end
    if showProgress; fprintf(1, '\nDone with isomap\n'); end
end

toc;

shiftRaw = dftresult(:, 3);
err = dftresult(:, 1);

%% Prepare output

try
    isomap.shift = reshape(shiftRaw, rs.navDims);
    isomap.err = reshape(err, rs.navDims);
catch
    % Lazy fix if the input is 1D rather than 2D
    isomap.shift = reshape(shiftRaw, [1 rs.navDims]);
    isomap.err = reshape(err, [1 rs.navDims]);
end
isomap.shiftRaw = isomap.shift;
isomap.interp = interp;
isomap.smoothData = smoothData;
isomap.diffOrder = diffOrder;
isomap.highPassWidth = highPassWidth;
% isomap.cropRight = abs(floor(min(isomap.shift(:))));
% isomap.cropLeft = abs(ceil(max(isomap.shift(:))));
% cropMask = true(size(ESI, 1), 1);
% cropMask((cropLeft+1) : (size(ESI, 1) - cropRight)) = false;

%% Waitbar function
function nUpdateWaitbar(~)
    waitbar(pp/rs.nPix, h);
    pp = pp + 1;
end

end

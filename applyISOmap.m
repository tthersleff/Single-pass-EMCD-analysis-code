function [ESIout, cropMask, shiftMask] = applyISOmap(ESI, isomap, varargin)
%APPLYISOMAP Applies an energy shift map to an ESI
%   applyISOmap takes an energy shift map (from the isomap variable) and
%   applies it to an entire spectrum image.  It uses sub-channel
%   interpolation through the shiftSubPx function.  Optionally, one can
%   'crop' the ESI to regions that were not shifted.  Alternatively, one
%   can retrieve the cropMask (length of the energy dimension of ESIout) or
%   the shiftMask (dimensions of the unfolded ESI) for further diagnostics.

%   (c) 2019 Thomas Thersleff, Stockholm University

%% Import the data

% Set default values
doCrop = false;
applyRaw = false;
invertShift = false; % Takes the negative of the shift map
useParallel = false;


% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'crop'
                doCrop = varargin{idx + 1};
            case 'interp'
                isomap.interp = varargin{idx + 1};
            case 'applyraw'
                applyRaw = varargin{idx + 1};
            case 'invertshift'
                invertShift = varargin{idx + 1};
            case 'parallel'
                useParallel = varargin{idx + 1};
        end
    end
end

if applyRaw
    isomap.shift = isomap.shiftRaw;
end

if invertShift
    isomap.shift = -isomap.shift;
end

% Matricize the ESI
[D, rs] = make2D(ESI, 1);
nE = rs.dims(rs.isSignal);

% Preallocate the new ESI
% expand = ceil(abs(max(isomap.shift(:)) - min(isomap.shift(:))));
expand = abs(max(round(isomap.shift(:))) - min(round(isomap.shift(:))));
nEout = nE + expand;
ESIout = zeros(nEout, rs.nPix);
SubPx = zeros(size(D));
shiftMask = true(size(ESIout));

% Get the shift values (fast shift)
shift = reshape(isomap.shift, rs.nPix, 1);
shift = (round(shift*10))/10;
shiftInt = round(shift);
shiftSub = shift - shiftInt;

iBegin = shiftInt - min(shiftInt) + 1;
iFinal = iBegin + nE - 1;

%% Now make the shift

if useParallel
    fprintf('\nParallel processing isn''t implemented yet. Sorry!\n');
    for idx = 1:rs.nPix
        SubPx(:, idx) = shift1D(D(:, idx), shiftSub(idx), true);
        ESIout(iBegin(idx):iFinal(idx), idx) = SubPx(:, idx);
        shiftMask(iBegin(idx):iFinal(idx), idx) = false;
    end
else
    for idx = 1:rs.nPix
        SubPx(:, idx) = shift1D(D(:, idx), shiftSub(idx), true);
        ESIout(iBegin(idx):iFinal(idx), idx) = SubPx(:, idx);
        shiftMask(iBegin(idx):iFinal(idx), idx) = false;
    end
end

ESIout = reshape(ESIout, [size(ESIout, 1), rs.navDims]);
cropMask = any(shiftMask, 2);

if doCrop
    ESIout(cropMask, :, :) = [];
end
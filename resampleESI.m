function [ resampled_ESI ] = resampleESI( ESI, new_dims, interpMethod )
%RESAMPLEESI Resamples an ESI to the new dimensions using
%griddedInterpolant and spline interpolation
%   resampleESI takes an input ESI and resamples it to the new dimensions
%   using Matlab's built-in griddedInterpolant function and applying a
%   spline interpolation.  It does not extrapolate.
%   IMPORTANT!  The ESI must not have any singleton dimensions.  The one
%   exception is if a 1d vector is passed, then I will account for this.
%   But if you did something like sum(ESI, 2) and get a vector with size 
%   [E 1 Y], then you need to pass: squeeze(ESI)
%   (c) 2019 Thomas Thersleff, Stockholm University

nDims = ndims(ESI);
dims = size(ESI);

%First, find out if "ESI" is a 1d vector.
%If the ESI has two dimensions and there is a singleton dimension:
if nDims == 2 && ~isempty(find(dims == 1))
    nDims = 1; %Set the number of dimensions to 1
end

F = griddedInterpolant(ESI);
if ~exist('interpMethod', 'var')
    F.Method = 'spline';
else
    F.Method = interpMethod;
end

    try 
        resize = new_dims ./ dims;
    catch ME
        if (strcmp(ME.identifier,'MATLAB:dimagree'))
            msg = sprintf('The new_dims input variable has %d dimensions, whereas your input ESI has %d dimensions',...
                ndims(new_dims), nDims);
            causeException = MException('MATLAB:myCode:dimensions',msg);
            ME = addCause(ME, causeException);
        end
        rethrow(ME)
    end

if nDims == 1
    v = linspace(1, max(dims), max(new_dims));
else
    for ind = 1:nDims
        v{ind, :} = linspace(1, dims(ind), round(dims(ind) .* resize(ind)));
    end
end

resampled_ESI = F(v);

end


function [ESIout, spikeESI] = applyXRaySpikeModel(ESI, spike)
%APPLYXRAYSPIKEMODEL Remove x-ray spikes from an EELS hyperspectral image
%   applyXRaySpikeModel simply applies the x-ray spike model (spike)
%   that is generated from the function getXRaySpikeModel.m

%   If tensorlab is installed, its sparse matrix representation is used.
%   Otherwise, Matlab's built-in sparse matrix representation is used.

%   Tensorlab can be downloaded from:
%   https://www.tensorlab.net/

% tensorlab; % Startup code for tensorlab package

try
    % Tensorlab implementation
    spikeESI = ful(spike);
catch
    % Matlab implementation
    spikeESI = full(spike);
end

ESIout = ESI - spikeESI;


end


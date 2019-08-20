function [ Vout ] = imGamma( Vin, gamma )
%IMGAMMA Apply a gamma curve to a colormap

A = 1;
Vout = A.*Vin.^(gamma);

end


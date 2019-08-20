function [y, shiftMask] = shift1D(S, shift, crop)
%SHIFT1D Shifts a 1D signal by a non-integer amount
%   shift1D shifts a 1D signal by an arbitrary amount.  It uses the Fourier
%   shift theorem to accomplish this, allowing for sub-channel accuracy.
%   It can account for both positive and negative shifts.  The maximum
%   value of shift is one half of the length of the signal vector.  It
%   avoids artifacts by mirroring the original spectrum.  By default, it
%   returns a vector of the same dimension as the input vector with all of
%   the shifted channels zeroed out.  It can alternatively return a vector
%   of length 2x the original length as well as a logical vector denoting
%   the channels that should be removed.
%   (c) 2019 Thomas Thersleff, Stockholm University


% Force column vector
S = S(:);

% Length of original vector
nE = length(S);
% mid = nE/2;

% Get shift rounded away from zero
delta = ceil(abs(shift)) * sign(shift);

% Make the zero vector. True means that the values should be zeroed out
z = [ones(nE, 1); zeros(nE, 1)];

% Mirror the vector
S = [flipud(S); S];
Nf = length(S);

% Function handle for the shift
fShift = @(sft) exp(-1i * 2 * pi * sft * [0:floor(Nf/2)-1 floor(-Nf/2):-1]' / Nf);

% Compute the linear phase shift
% W = exp(-1i * 2 * pi * shift * [0:floor(Nf/2)-1 floor(-Nf/2):-1]' / Nf); 
% Wz = exp(-1i * 2 * pi * delta * [0:floor(Nf/2)-1 floor(-Nf/2):-1]' / Nf); 
W = fShift(shift);
Wz = fShift(delta);

% Force conjugate symmetry.  Always needed b/c mirroring forces the vector
% to be even in length
W(Nf/2+1) = real(W(Nf/2+1));

% Compute the phase shift
X = fft(S);
Y = X .* W;
Z = fft(z) .* Wz;

% Invert the FFT.
y = ifft(Y);
zf = ifft(Z);

% Ensure that the left half of the shift vector is removed
% z(1:nE) = 1;
% z = logical(round(z));

% There should be no imaginary component (for real input
% signals) but due to numerical effects some remnants remain.
if isreal(S)
    y = real(y);
end
try
    zf = logical(round(zf));
catch
    zf = zf;
end
zf(1:nE) = true;
y(zf) = 0;
    
if crop
    y(1:nE) = [];
else
    Wf = fShift(-Nf/4);
    Yf = fft(y) .* Wf;
    y = real(ifft(Yf));
    shiftMask = real(ifft(fft(zf) .* Wf));
    shiftMask = logical(round(shiftMask));
end


end


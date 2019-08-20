function [sse, mLmS, out] = objMLMS(e, plus, minus, f, win, varargin)
%OBJMLMS The objective function for mLmS optimization
%   Requires prepEMCDoptim to be run first

%% Set default values
preModel = 'Power Law'; % 'Power Law' or 'LCPL'
postModel = 'Integral'; % 'Integral' or 'Linear' or 'Power'
SScalculation = 'Coupled'; % 'Coupled' or 'Decoupled' or 'Cumulative'
sigModel = 'Gaussian'; % 'Gaussian' or 'Skewed' or 'GLmix'
forceOppositeSign = true;
% doDeblur = false;
% doShift = false;

%% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'premodel'
                preModel = varargin{idx + 1};
            case 'postmodel'
                postModel = varargin{idx + 1};
            case 'sscalculation'
                SScalculation = varargin{idx + 1};
            case 'sigmodel'
                sigModel = varargin{idx + 1};
            case 'forceoppositesign'
                forceOppositeSign = varargin{idx + 1};
%             case 'deblur'
%                 doDeblur = varargin{idx + 1};
%             case 'shift'
%                 doShift = varargin{idx + 1};
        end
    end
end

%% Constrain peaks to be of opposite sign if requested

if forceOppositeSign
    if f(13)*f(19) > 0 % If both have the same sign
        f(19) = -f(19); % Flip the sign on L2
    end
end

%% Model and remove background

if strcmpi(preModel, 'Power Law')
    BGplus = f(1) .* e .^ f(2);
    BGminus = f(3) .* e .^ f(4);
elseif strcmpi(preModel, 'LCPL')
    BGplus = (f(1:2) * [e(:).^-3 e(:).^-2]')';
    BGminus = (f(3:4) * [e(:).^-3 e(:).^-2]')';
else
    error('Unrecognized pre-edge model');
end

preplus = plus - BGplus;
preminus = minus - BGminus;

%% Perform deblurring and shifting
preminus = preminus - f(25).*gradient(gradient(smooth(preminus)));
preminus = shift1D(preminus, f(26), true);

%% Post-edge normalization

if strcmpi(postModel, 'Integral')
    preminus = preminus .* f(9);
elseif strcmpi(postModel, 'Linear')
    ln = f(9) + e .* f(10);
    preminus = preminus .* ln;
elseif strcmpi(postModel, 'Power')
    ln = exp(f(9)) .* e .^ f(10);
    preminus = preminus .* ln;
else
    error('Unrecognized post-edge model');
end

%% Normalize to max value of L3

mp = max(preplus(makeWin(e, [707 712])));
mm = max(preminus(makeWin(e, [707 712])));
mx = max(mp, mm);
preplus = preplus ./ mx;
preminus = preminus ./ mx;

%% Get EMCD

EMCD = preplus - preminus;

%% Extract EMCD fit from fitting parameters

if strcmpi(sigModel, 'Gaussian')
    L3 = @(x) x(13).*exp(-((e-x(14))./x(15)).^2);
    L2 = @(x) x(19).*exp(-((e-x(20))./x(21)).^2);
elseif strcmpi(sigModel, 'Skewed')
    L3 = @(x) x(13).*exp((x(14)/2).*(2*x(15) + x(14)*x(16)^2 - 2.*e)).*erfc((x(15)+x(14)*x(16)^2-e)/sqrt(2)*x(16));
    L2 = @(x) x(19).*exp((x(20)/2).*(2*x(21) + x(20)*x(22)^2 - 2.*e)).*erfc((x(21)+x(20)*x(22)^2-e)/sqrt(2)*x(22));
elseif strcmpi(sigModel, 'GLmix')
    L3 = @(x) x(13).*x(17)*(1+((e-x(16))./(x(14)+x(15).*(e-x(16)))).^2).^-1 + ...
              x(13).*(1-x(17)).*exp(log(0.5).*(((e-x(16))./(x(14)+x(15).*(e-x(16)))).^2));
    L2 = @(x) x(19).*x(17).*(1+((e-x(22))./(x(20)+x(21).*(e-x(22)))).^2).^-1 + ...
              x(19).*(1-x(17)).*exp(log(0.5).*(((e-x(22))./(x(20)+x(21).*(e-x(22)))).^2));
else
    error('Unrecognized EMCD signal model');
end

EMCDsig = L3(f) + L2(f);

%% Calculate mLmS
fK = @(p, q) (2.*q)./(9.*p - 6.*q);
% Removing trapz to try to increase speed
p = trapz(e, L3(f));
q = trapz(e, L2(f)) + p;
% p = sum(L3(f));
% q = sum(L2(f)) + p;

mLmS = fK(p, q);

%% Calculate the squared sum of errors (SSE)

% SSE regions
preBegin = find(e >= win.pre(1), 1);
preEnd = find(e <= win.pre(2), 1, 'last');
% postBegin = find(e > win.pre(2), 1, 'first');
postBegin = find(e > win.post(1), 1, 'first');
postEnd = find(e < win.post(2), 1, 'last');

L3Begin = find(e > win.pre(2), 1, 'first');
L3End = find(e < 718, 1, 'last');
L2Begin = find(e > 718, 1, 'first');
L2End = find(e < win.post(1), 1, 'last');

% Residuals
res = EMCD - EMCDsig;

if strcmpi(SScalculation, 'Coupled')
    sse = sum(res(preBegin:postEnd) .^ 2);
elseif strcmpi(SScalculation, 'Decoupled')
    resPlus = preplus(preBegin:preEnd);
    resMinus = preminus(preBegin:preEnd);
    
    sse = ...
      10 .* sum(res(postBegin:postEnd) .^ 2) + ...
      1 .* sum(res(L3Begin:L3End) .^ 2) + ...
      1 .* sum(res(L2Begin:L2End) .^ 2) + ...
      1 .* sum(resPlus .^ 2) + ...
      1 .* sum(resMinus .^ 2) + ...
      0 .* sum(EMCDsig(L2Begin:L2End).^2) + ...
      0 .* sum(EMCDsig(L3Begin:L3End).^2) + ...
      10 .* sum(EMCDsig(postBegin:postEnd).^2);
  %       1 .* sum(res(preBegin:preEnd) .^ 2);
elseif strcmpi(SScalculation, 'Cumulative')
    resPlus = preplus(preBegin:preEnd);
    resMinus = preminus(preBegin:preEnd);
    
    cumSig = cumsum(EMCD(preBegin:postEnd));
    cumFit = cumsum(EMCDsig(preBegin:postEnd));
    
    cumRes = cumSig - cumFit;
    
    sse = sum(cumRes .^ 2);
else
    error('Unrecognized SSE coupling mode');
end


%% Return all data if requested
if nargout == 3
    out.preplus = preplus;
    out.preminus = preminus;
    out.BGplus = BGplus;
    out.BGminus = BGminus;
    out.mp = mp;
    out.mm = mm;
    out.mx = mx;
    out.EMCD = EMCD;
    out.EMCDsig = EMCDsig;
    out.p = p;
    out.q = q;
    out.mLmS = mLmS;
    out.e = e;
    out.f = f;
    out.plus = plus;
    out.minus = minus;
    out.res = res;
    out.sse = sse;
    out.std = std(res(preBegin:postEnd));
    out.SNRL3 = 10*log10(norm(EMCDsig(L3Begin:L3End)) ./ norm(res(L3Begin:L3End)));
    out.SNRL2 = 10*log10(norm(EMCDsig(L2Begin:L2End)) ./ norm(res(L2Begin:L2End)));
end

end


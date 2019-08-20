function [F, b] = prepEMCDoptim(e, Plus, Minus, window, varargin)
%PREPEMCDOPT Prepares for the optimization function
%   F(1:8) --> Pre-edge BG parameters
%   F(9:12) --> Post-edge normalization parameters
%   F(13:18) --> EMCD L3 signal parameters
%   F(19:24) --> EMCD L2 signal parameters

%% Set default values
preModel = 'Power Law'; % 'Power Law' or 'LCPL'
postModel = 'Integral'; % 'Integral' or 'Linear' or 'Power'
SScalculation = 'Coupled'; % 'Coupled' or 'Decoupled' or 'Cumulative'
sigModel = 'Gaussian'; % 'Gaussian' or 'Skewed' or 'GLmix'
doDeblur = false;
doShift = false;

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
            case 'deblur'
                doDeblur = varargin{idx + 1};
            case 'shift'
                doShift = varargin{idx + 1};
        end
    end
end

% Make the pre and post masks
% pre = makeWin(e, window.pre);
post = makeWin(e, window.post);

%% Step 1: get the pre-edge BG

if strcmpi(preModel, 'Power Law')
    % Parameters
    [Pfit, preplus] = afit(e, Plus, window.pre);
    [Mfit, preminus] = afit(e, Minus, window.pre);
    F(1) = Pfit(1);
    F(2) = Pfit(2);
    F(3) = Mfit(1);
    F(4) = Mfit(2);
    F(5:8) = 0;
    % Bounds
    if strcmpi(SScalculation, 'Coupled')
        Arange = 0.00;
        Rrange = 0.00;
        % Lower bounds
        b.lb(1) = F(1) - abs((F(1) * Arange));
        b.lb(2) = F(2) - abs((F(2) * Rrange));
        b.lb(3) = F(3) - abs((F(3) * Arange));
        b.lb(4) = F(4) - abs((F(4) * Rrange));    
        % Upper bounds     
        b.ub(1) = F(1) + abs((F(1) * Arange));
        b.ub(2) = F(2) + abs((F(2) * Rrange));
        b.ub(3) = F(3) + abs((F(3) * Arange));
        b.ub(4) = F(4) + abs((F(4) * Rrange));
    elseif strcmpi(SScalculation, 'Decoupled')
%         Arange = 1;
%         Rrange = 1;
%         b.lb(1:4) = -Inf;
        b.lb([1 3]) = 0;
        b.lb([2 4]) = -5;
%         b.ub(1:4) = +Inf;
        b.ub([1 3]) = +Inf;
        b.ub([2 4]) = 0;
    elseif strcmpi(SScalculation, 'Cumulative')
        b.lb(1:4) = -Inf;
        b.ub(1:4) = +Inf;
    else
        error('Unrecognized SSE coupling mode');
    end

    b.lb(5:8) = 0;
    b.ub(5:8) = 0;
    
elseif strcmpi(preModel, 'LCPL')
    [Pfit, ~, preplus] = LCPL(e, Plus, window.pre);
    [Mfit, ~, preminus] = LCPL(e, Minus, window.pre);
    F(1) = Pfit.c(1);
    F(2) = Pfit.c(2);
    F(3) = Mfit.c(1);
    F(4) = Mfit.c(2);
    F(5:8) = 0;
    % Bounds
    if strcmpi(SScalculation, 'Coupled')
        Arange = 0.00;
        Rrange = 0.00;
    elseif strcmpi(SScalculation, 'Decoupled')
        Arange = 1000;
        Rrange = 1000;
    elseif strcmpi(SScalculation, 'Cumulative')
        Arange = 1000;
        Rrange = 1000;
    else
        error('Unrecognized SSE coupling mode');
    end
    % Lower bounds
    b.lb(1) = F(1) - abs((F(1) * Arange));
    b.lb(2) = F(2) - abs((F(2) * Rrange));
    b.lb(3) = F(3) - abs((F(3) * Arange));
    b.lb(4) = F(4) - abs((F(4) * Rrange));
    b.lb(5:8) = 0;
    % Upper bounds
    b.ub(1) = F(1) + abs((F(1) * Arange));
    b.ub(2) = F(2) + abs((F(2) * Rrange));
    b.ub(3) = F(3) + abs((F(3) * Arange));
    b.ub(4) = F(4) + abs((F(4) * Rrange));
    b.ub(5:8) = 0;
else
    error('Unrecognized pre-edge model');
end

%% Step 2: Get the post-edge normalization

if strcmpi(postModel, 'Integral')
    top = sum(preplus(post));
    bottom = sum(preminus(post));
    N = top/bottom;
    F(9) = abs(N); %Needed for very noisy spectra, but biases result
    F(10:12) = 0;
    % Bounds
    b.lb(9) = F(9) - abs(100*F(9));
    b.ub(9) = F(9) + abs(100*F(9));
    b.lb(10:12) = 0;
    b.ub(10:12) = 0;
elseif strcmpi(postModel, 'Linear')
    D = preplus(post) ./ preminus(post);
    N = robustfit(e(post), D);
    F(9) = N(1);
    F(10) = N(2);
    F(11:12) = 0;
    % Bounds
    b.lb(9) = F(9) - abs(100*F(9));
    b.ub(9) = F(9) + abs(100*F(9));
    b.lb(10) = F(10) - abs(100*F(10));
    b.ub(10) = F(10) + abs(100*F(10));
    b.lb(11:12) = 0;
    b.ub(11:12) = 0;
elseif strcmpi(postModel, 'Power')
    D = preplus(post) ./ preminus(post);
    N = robustfit(log(e(post)), log(D));
    F(9) = N(1);
    F(10) = N(2);
    F(11:12) = 0;
    % Bounds
    b.lb(9) = F(9) - abs(100*F(9));
    b.ub(9) = F(9) + abs(100*F(9));
    b.lb(10) = F(10) - abs(100*F(10));
    b.ub(10) = F(10) + abs(100*F(10));
    b.lb(11:12) = 0;
    b.ub(11:12) = 0;
else
    error('Unrecognized post-edge model');
end

%% Step 3: Get the EMCD signal model

if strcmpi(sigModel, 'Gaussian')
    % L3
    F(13) = 0.05;
    F(14) = 709.8;
    F(15) = 2;
    F(16:18) = 0; % Unused
    % L2
    F(19) = -0.05;
    F(20) = 722.1;
    F(21) = 3;
    F(22:24) = 0; % Unused
    
    % Bounds
    % L3
    b.lb(13) = -1; 
    b.ub(13) = 1;
    b.lb(14) = 706;
    b.ub(14) = 715;
    b.lb(15) = 0.5;
    b.ub(15) = 2.0;
    b.lb(16:18) = 0;
    b.ub(16:18) = 0;
    % L2
    b.lb(19) = -1;
    b.ub(19) = 1;
    b.lb(20) = 718;
    b.ub(20) = 726;
    b.lb(21) = 1.0;
    b.ub(21) = 5.0;
    b.lb(22:24) = 0;
    b.ub(22:24) = 0;
    
elseif strcmpi(sigModel, 'Skewed')
    % L3
    F(13) = 0.1;
    F(14) = 1.5;
    F(15) = 710;
    F(16) = 1.3;
    F(17:18) = 0; % Unused
    % L2
    F(19) = -0.05;
    F(20) = 0.5;
    F(21) = 722;
    F(22) = 0.7;
    F(23:24) = 0; % Unused
    
    % Bounds
    % L3
    b.lb(13) = -1; 
    b.ub(13) = 1;
    b.lb(14) = 0.5;
    b.ub(14) = 2.0;
    b.lb(15) = 700;
    b.ub(15) = 712;
    b.lb(16) = 0.5;
    b.ub(16) = 2.0;
    b.lb(17:18) = 0;
    b.ub(17:18) = 0;
    % L2
    b.lb(19) = -1; 
    b.ub(19) = 1;
    b.lb(20) = 0.1;
    b.ub(20) = 2.0;
    b.lb(21) = 716;
    b.ub(21) = 725;
    b.lb(22) = 0.1;
    b.ub(22) = 2.0;
    b.lb(23:24) = 0;
    b.ub(23:24) = 0;
elseif strcmpi(sigModel, 'GLmix')
    % L3
    F(13) = 0.05;
    F(14) = 1.3;
    F(15) = 0.1;
    F(16) = 710;
    F(17) = 0.33; %Lmix
    F(18) = 0; % Unused
    % L2
    F(19) = -0.03;
    F(20) = 2.3;
    F(21) = 0.17;
    F(22) = 723;
    F(23:24) = 0; % Unused
    
    % Bounds
    % L3
    b.lb(13) = -1; 
    b.ub(13) = 1;
    b.lb(14) = 0; % 0
    b.ub(14) = 5.0;
    b.lb(15) = 0;
    b.ub(15) = 5.0;
    b.lb(16) = 706;
    b.ub(16) = 712;
    b.lb(17) = 0; % Lmix
    b.ub(17) = 1; % Lmix
    b.lb(18) = 0; % Unused
    b.ub(18) = 0; % Unused
    % L2
    b.lb(19) = -1; 
    b.ub(19) = 1;
    b.lb(20) = 0;
    b.ub(20) = 5.0; %5.0
    b.lb(21) = 0.14;
    b.ub(21) = 5.0; 
    b.lb(22) = 718;
    b.ub(22) = 725;
    b.lb(23:24) = 0;
    b.ub(23:24) = 0;
else
    error('Unrecognized EMCD signal model');
end

if doDeblur
    F(25) = 0;
    b.lb(25) = -4;
    b.ub(25) = +4;
else
    F(25) = 0;
    b.lb(25) = 0;
    b.ub(25) = 0;
end

if doShift
    F(26) = 0;
    b.lb(26) = -3;
    b.ub(26) = +3;
else
    F(26) = 0;
    b.lb(26) = 0;
    b.ub(26) = 0;
end

end


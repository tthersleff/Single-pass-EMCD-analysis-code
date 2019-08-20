function [fo, mLmS, fStruct, film] = optimMLMS(e, plus, minus, win, varargin)
%OPTIMMLMS Optimization wrapper for mLmS modelling
%   optimMLMS is an optimization routine intended to improve the
%   reliability of EMCD signal extraction from two chiral signals.  It is
%   detailed in the paper by Thersleff et al. and, if it is used, should be
%   cited.  A DOI with more information is here:
%   10.5281/zenodo.3361582

%   (c) 2019 Thomas Thersleff, Stockholm University

%% Set default values
preModel = 'Power Law'; % 'Power Law' or 'LCPL'
postModel = 'Integral'; % 'Integral' or 'Linear' or 'Power'
SScalculation = 'Coupled'; % 'Coupled' or 'Decoupled' or 'Cumulative'
sigModel = 'Gaussian'; % 'Gaussian' or 'Skewed'
forceOppositeSign = true;
algorithm = 'interior-point'; % 'interior-point' 'active-set' 'sqp'
pauseTime = 0.005;
xlim = [660 860];
tol = 1e-12;
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
            case 'forceoppositesign'
                forceOppositeSign = varargin{idx + 1};
            case 'pausetime'
                pauseTime = varargin{idx + 1};
            case 'algorithm'
                algorithm = varargin{idx + 1};
            case 'xlim'
                xlim = varargin{idx + 1};
            case 'functiontolerance'
                tol = varargin{idx + 1};
            case 'deblur'
                doDeblur = varargin{idx + 1};
            case 'shift'
                doShift = varargin{idx + 1};
        end
    end
end

film = [];

%% Step 1: Get starting conditions

[fi, b] = prepEMCDoptim(e, plus, minus, win, ...
    'preModel', preModel, ...
    'postModel', postModel, ...
    'SScalculation', SScalculation, ...
    'sigModel', sigModel, ...
    'deblur', doDeblur, ...
    'shift', doShift);

%% Initialize graph if requested

if nargout == 4
    [~, mLmS, fStructInit] = objMLMS(e, plus, minus, fi, win, ...
        'preModel', preModel, ...
        'postModel', postModel, ...
        'SScalculation', SScalculation, ...
        'sigModel', sigModel, ...
        'forceOppositeSign', forceOppositeSign);
%     h = plotEMCDfit(fStructInit, 'xlim', xlim);
    h = plotEMCDfitCum(fStructInit, 'xlim', xlim);
%     title(h.ax(3), sprintf('m_L/m_S = %0.4f', mLmSinit));
%     title(h.ax(2), sprintf('m_L/m_S = %0.4f', fK(p, q)));

    film = [film getframe(h.fig)];
    opts = optimoptions(@fmincon,'OutputFcn',@outfun,... 
        'Display','iter',...
        'Algorithm', algorithm);
    
else
    opts = optimoptions(@fmincon, ... 
        'Display','none',...
        'Algorithm', algorithm);
end

%% Final changes to optimset

% opts.ConstraintTolerance = 1e-12;
% opts.OptimalityTolerance = 1e-12;
% opts.FunctionTolerance = 1e-12; % Good choice
% opts.FunctionTolerance = 1e-8; % Faster convergence
opts.FunctionTolerance = tol;
opts.UseParallel = false;

%% Set the objective function

fun = @(ff) objMLMS(e, plus, minus, ff, win, ...
        'preModel', preModel, ...
        'postModel', postModel, ...
        'SScalculation', SScalculation, ...
        'sigModel', sigModel, ...
        'forceOppositeSign', forceOppositeSign);

%% Optimization
% [fo, ~, ~, ~, ~, ~, hessian] = fmincon(fun, fi, [], [], [], [], b.lb, b.ub, [], opts);
% try
    fo = fmincon(fun, fi, [], [], [], [], b.lb, b.ub, [], opts);
% catch
% %     Sometimes the optimization fails if the spectra are truly bad.  This
% %     simply sets the output parameters to the estimated ones.  This is
% %     useful when running this code in a loop over very noisy data where
% %     EMCD may not be expected to be present.
%     fo = fi;
% end


%% Plotting function (if requested)

    function stop = outfun(x, ~, state)
         stop = false;
 
         switch state
             case 'init'
                % Code
                pause(pauseTime);
             case 'iter'
                [~, mLmS, fStruct] = objMLMS(e, plus, minus, x, win, ...
                    'preModel', preModel, ...
                    'postModel', postModel, ...
                    'SScalculation', SScalculation, ...
                    'sigModel', sigModel, ...
                    'forceOppositeSign', forceOppositeSign);
                Lobj = h.pl;
                set(findobj(Lobj, 'DisplayName', 'BG Plus'), 'YData', fStruct.BGplus);
                set(findobj(Lobj, 'DisplayName', 'BG Minus'), 'YData', fStruct.BGminus);
                set(findobj(Lobj, 'DisplayName', 'Plus BGsub'), 'YData', fStruct.preplus);
                set(findobj(Lobj, 'DisplayName', 'Minus BGsub'), 'YData', fStruct.preminus);
                set(findobj(Lobj, 'DisplayName', 'EMCD signal'), 'YData', fStruct.EMCD);
                set(findobj(Lobj, 'DisplayName', 'EMCD fit'), 'YData', fStruct.EMCDsig);
                set(findobj(Lobj, 'DisplayName', 'EMCD cumsum'), ...
                    'YData', cumsum(fStruct.EMCD(makeWin(e, xlim))));
                set(findobj(Lobj, 'DisplayName', 'EMCD fit cumsum'), ...
                    'YData', cumsum(fStruct.EMCDsig(makeWin(e, xlim))));
                title(findobj(h.ax, 'Tag', 'Cumsum Plot'), sprintf('m_L/m_S = %0.4f', mLmS));
                title(findobj(h.ax, 'Tag', 'Norm Plot'), sprintf('L_3 = %2.2f %%, L_2 = %2.2f %%', abs(fStruct.f(13)*100), abs(fStruct.f(19)*100)));
                title(findobj(h.ax, 'Tag', 'EMCD Plot'), sprintf('SNR L_3 = %2.2f dB, L_2 = %2.2f dB', fStruct.SNRL3, fStruct.SNRL2));
                drawnow;
                film = [film getframe(h.fig)];
                pause(pauseTime);

             case 'done'
                 hold off
             otherwise
        end
    end 

%% Force opposite sign for output

if forceOppositeSign
    if fo(13)*fo(19) > 0 % If both have the same sign
        fo(19) = -fo(19); % Flip the sign on L2
    end
end


%% Get additional output if requested

if nargout >= 2
    [~, mLmS, fStruct] = objMLMS(e, plus, minus, fo, win, ...
                    'preModel', preModel, ...
                    'postModel', postModel, ...
                    'SScalculation', SScalculation, ...
                    'sigModel', sigModel, ...
                    'forceOppositeSign', forceOppositeSign);
end

end



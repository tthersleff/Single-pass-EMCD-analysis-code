function [ data, gui ] = viewESI( e, dimSignal, varargin )
%viewESI Visualization of hyperspectral EELS datacubes
%   GUI visualization of hyperspectral datacubes, intended for EELS and EDX
%   analysis.  
%   Requires the GUI Layout Toolbox:
%   http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
%   (c) Thomas Thersleff, 2016, Uppsala University
%
%
% INPUT:    e           [1 x L]     --> Energy vector for the EELS datacube.  
%                                       Typically calibrated in units of eV
%           varargin    [M x N x L] --> The hyperspectral datacubes.  The
%                                       energy axis needs to be dimension 3 
%                                       (referred to here as "L")  There is
%                                       no fundamental limit to how many of
%                                       these datacubes can be passed.
%                                       They will all be shown
%                                       simultaneously on the graph
%
% OUTPUT:   data            The data struct
%           data.dims3D     Dimensions of the 3D datacube inputs
%           data.dims4D     The datacubes are written to a 4D array to simplify
%                           their treatment.  Dimensions are [A x M x N x L]
%                           where "A" is now the index of which ESI is being
%                           used (essentially the order in which the ESI are
%                           passed to the function).
%           data.ESI        The 4D ESI datacube used to contain all ESI
%                           datacubes passed to the function
%           data.ROI        The ROI for the summation
%           data.map        The map used to define the spatial coordinates.  At
%                           some point this will be updated so that you can
%                           pass a different image (such as an HAADF).  The
%                           default is to simply take the sum of the lowest
%                           order ESI datacube along the L dimension (energy
%                           loss)
%           data.S          The Spectra to be plotted.  Dimensions [A x L]
%
%           gui             The figure struct.  Used to build the figure but
%                           otherwise not so important.
%
%
% ACRONYMS: ESI     Electron Spectrum Image
%           SI      Spectrum Image
%           EELS    Electron Energy Loss Spectroscopy
%           EDX     Energy Dispersive X-Ray Spectroscopy.  Alternatively
%                   known as EDS and EDXS.
%           ROI     Region Of Interest.  Denotes the location on the map
%                   that is used for taking the spectral mean that is to be
%                   displayed in the graph.  Defined using Matlab's
%                   built-in imrect function.

data = initializeData();


%% Build Initial Data Struct

    function data = initializeData()
        % Find which dimension is the signal
        if ~exist('dimSignal', 'var')
            dimSignal = 1;
        end
        isSignal = false(1, ndims(varargin{1}));
        isSignal(dimSignal) = true;
        allDims = 1:ndims(varargin{1});
        navDims = allDims(~isSignal);
        permOrd = [navDims dimSignal];
        
        % Get all the ESI datacubes and write them to a single 4D matrix, data.ESI       
        if isempty(varargin)
            error('myfuns:viewESI:TooManyInputs', ...
                'This function requires at least 1 ESI after the energy vector');
        end

        % Preallocation for speed
        data.dims3D = size(varargin{1});
        data.dims3D = data.dims3D(permOrd);
        data.dims4D = [length(varargin) data.dims3D(1) data.dims3D(2) data.dims3D(3)];

        data.ESI = zeros(data.dims4D);

        % Get all ESI datacubes
        for ind = 1:length(varargin)
            data.ESI(ind,:,:,:) = permute(varargin{ind}, permOrd);
        end

        data.ROI = [round(data.dims4D(3)/2) round(data.dims4D(2)/2) 1 1]; %place sample rectangle in the middle
        %NOTE: The ROI is placed using x as the first number and y as the second.
        %This is different from matrix indexing, which uses y as the first number
        %(column) and x as the second.


        data.map = squeeze(sum(data.ESI(1,:,:,:),4));
        
    end % initializeData

%% Build the figure

% Get the correct figure dimensions (so complicated!)
gui.screen = get(0,'screensize'); %How big is the screen?
gui.figuresize = gui.screen;
gui.figuresize(3) = 0.75*gui.screen(3); % Set figure width
gui.figuresize(4) = 0.75*gui.screen(4); % Set figure height.
gui.fig = figure('Color',[1 1 1],...
                 'OuterPosition',gui.figuresize,...
                 'Units','normalized',...
                 'Visible','off',...
                 'NumberTitle', 'off', ...
                 'Name','View ESI'); % Finally make the figure

movegui(gui.fig,'center');
set(gui.fig,'Visible','on'); %Make the figure visible after moving

% If you prefer a different color order, define it here
gui.colorOrder = [0.5 0.5 0.5;...
                  0.635,  0.078,  0.184;...
                  0.466,  0.674,  0.188;...
                  0.929,  0.694,  0.125;...
                  0.301,  0.745,  0.933];

% Make the main layout
gui.mainLayout = uix.HBoxFlex('Parent', gui.fig, 'Spacing', 3);
gui.panelMaps = uix.BoxPanel('Parent', gui.mainLayout, 'Title', 'Maps');

% Layout for the graphs panel
gui.layoutGraphs = uix.VBoxFlex('Parent', gui.mainLayout, 'Spacing', 3);
gui.panelMainGraph = uix.BoxPanel('Parent', gui.layoutGraphs, 'Title', 'Spectra');
gui.panelFormat = uix.Panel('Parent', gui.layoutGraphs);
gui.panelFormatSplit = uix.HBoxFlex('Parent',gui.panelFormat, 'Spacing', 3);

% Layout for the format axis panel
axisNames = {'X'; 'Y'};
    for idx = 1:2
        gui.panelFormatAxis(idx) = uix.BoxPanel('Parent',gui.panelFormatSplit,...
                                              'Title', strcat({'Format' },axisNames{idx}));
        gui.gridFormatAxis(idx) = uix.GridFlex('Parent',gui.panelFormatAxis(idx),...
                                               'DividerMarkings','off');

        gui.grids{idx,1} = uicontrol('Style', 'text',...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'String', strcat(axisNames{idx}, ' Lim Lower'));   
                                 
        gui.grids{idx,2} = uicontrol('Style', 'edit',...
                                     'Tag', strcat('Scale',axisNames{idx},'Lower'),...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'Callback',@onFormatAxis);

        gui.grids{idx,3} = uicontrol('Style', 'text',...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'String', strcat(axisNames{idx}, ' Lim Upper')); 
                                 
        gui.grids{idx,4} = uicontrol('Style', 'edit',...
                                     'Tag', strcat('Scale',axisNames{idx},'Upper'),...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'Callback',@onFormatAxis);

        gui.grids{idx,5} = uicontrol('Style', 'text',...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'String', strcat(axisNames{idx}, ' Log')); 

        gui.grids{idx,6} = uicontrol('Style', 'checkbox',...
                                     'Tag', strcat('Log',axisNames{idx}),...
                                     'Parent', gui.gridFormatAxis(idx),...
                                     'Callback',@onFormatAxis);
              
        set(gui.gridFormatAxis(idx), 'Widths',[-1 -1 -1], 'Heights', [-1 -1]);
              
    end
                         

set(gui.mainLayout, 'Widths',[-1 -2]);
set(gui.layoutGraphs, 'Heights',[-10 -1]);
% set(gui.layoutFormatAxis, 'Widths', [-1 -1 -1], 'MinimumWidths',[100 100 100]);

% Draw the axes
gui.axMap = axes('Parent', gui.panelMaps,...
    'ActivePositionProperty', 'Position',...
    'xtick',[],'ytick',[],...
    'YDir','reverse'); %YDir reverse is apparently needed to keep the direction correct

    %'OuterPosition',[0 .25 0.3 .75],...
    
gui.axSpectra = axes('Parent', uicontainer('Parent', gui.panelMainGraph),...
    'colororder', gui.colorOrder,...
    'ActivePositionProperty', 'outerposition',...
    'fontsize',12,...
    'xgrid','on','ygrid','on');

    %'OuterPosition',[0.3 .25 0.7 .75],...
    

    
colormap(gui.axMap, 'gray');
hold(gui.axMap,'on');
axis(gui.axMap,'off');



xlabel(gui.axSpectra, 'Energy [eV]'); ylabel(gui.axSpectra, 'Intensity [counts]');
hold(gui.axSpectra,'on');


% Plot the initial data
axes(gui.axMap);
gui.dataMap = imagesc(data.map); axis image;


%axes(gui.axSpectra)
gui.dataSpectrum = plot(gui.axSpectra, e, getSpectraFromROI(data.ESI, data.ROI));
gui.Lims(1,:) = get(gui.axSpectra, 'xlim'); gui.Lims(2,:) = get(gui.axSpectra, 'ylim');
set(gui.grids{1,2}, 'String', num2str(gui.Lims(1,1)));
set(gui.grids{1,4}, 'String', num2str(gui.Lims(1,2)));
set(gui.grids{2,2}, 'String', num2str(gui.Lims(2,1)));
set(gui.grids{2,4}, 'String', num2str(gui.Lims(2,2)));
% set(gui.formatAxis.xUp, 'String', num2str(gui.Lims(1,2)));



% Set the ROI rectangle

gui.imrect = imrect(gui.axMap, data.ROI);
setResizable(gui.imrect,1);
setColor(gui.imrect,[0 0 1]);

%Add callback to the ROI position
addNewPositionCallback(gui.imrect,@(p) updateROI(p));

% Constrain the ROI position to the axis limits
gui.constrainPosition = makeConstrainToRectFcn('imrect',...
    get(gui.axMap,'XLim'),get(gui.axMap,'YLim'));
setPositionConstraintFcn(gui.imrect, gui.constrainPosition); 

updateROI(data.ROI);



%% Callback functions

    function updateROI(p)
        data.ROI = round(p);
        title(gui.axMap, mat2str(data.ROI,3));
        data.S = getSpectraFromROI(data.ESI, data.ROI);
        updatePlots();
    end

    function updatePlots()
        for ind = 1:data.dims4D(1)
            set(gui.dataSpectrum(ind), 'ydata', data.S(ind,:));
        end
    end

    function onFormatAxis(source,~)
        Obj2Change = get(source,'Tag');
        
        switch Obj2Change
            
            case 'LogX'
                x = get(source, 'Value');
                    if x
                        gui.axisScale{1} = 'log';
                    else
                        gui.axisScale{1} = 'linear';
                    end
                set(gui.axSpectra,'xscale', gui.axisScale{1});
                
            case 'ScaleXLower'
                
                x = get(source, 'String');
                gui.Lims(1,1) = str2double(x);
                set(gui.axSpectra,'xlim',gui.Lims(1,:));
                
            case 'ScaleXUpper'
                
                x = get(source, 'String');
                gui.Lims(1,2) = str2double(x);
                set(gui.axSpectra,'xlim',gui.Lims(1,:));
                
            case 'LogY'
                x = get(source, 'Value');
                    if x
                        gui.axisScale{2} = 'log';
                    else
                        gui.axisScale{2} = 'linear';
                    end
                set(gui.axSpectra,'yscale', gui.axisScale{2});
                
            case 'ScaleYLower'
                
                x = get(source, 'String');
                gui.Lims(2,1) = str2double(x);
                set(gui.axSpectra,'ylim',gui.Lims(2,:));
                
            case 'ScaleYUpper'
                
                x = get(source, 'String');
                gui.Lims(2,2) = str2double(x);
                set(gui.axSpectra,'ylim',gui.Lims(2,:));
                
        end
        
    end

%% Other functions

    function S = getSpectraFromROI(ESI, ROI)

          S = squeeze2(mean(mean(ESI(:,...
                  ROI(2):ROI(2) + ROI(4) - 1,...
                  ROI(1):ROI(1) + ROI(3) - 1,...
                  :),3),2));
    end

    function b = squeeze2(a)
    % Found on the Mathworks site, doesn't return the transpose of the squeezed
    % matrix.  
    % http://se.mathworks.com/matlabcentral/answers/96861-why-does-the-squeeze-function-return-the-transpose-of-the-expected-output-in-matlab-7-2-r2006a
        S = size(a);
        if S(1) == 1
            b = squeeze(a)';
        else
            b = squeeze(a);
        end
    end
        
end


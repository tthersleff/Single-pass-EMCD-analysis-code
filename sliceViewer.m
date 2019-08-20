function h = sliceViewer( ESI, doPermutation )
%SLICEVIEWER Use a scroll bar to view all slices of a 3D matrix
%   sliceViewer is a simple script to allow you to scroll through the
%   images in a 3D matrix.  The first dimension must be the index
%   dimension, and the final two dimensions are then the spatial
%   dimensions.  The function imagesc is used to consistently scale the
%   image, and a color bar is added
%   (c) 2018 Thomas Thersleff, Stockholm University

if ~exist('doPermutation', 'var')
    doPermutation = 0;
end

if doPermutation
    ESI = permute(ESI, [3 1 2]);
end

%% Build the figure

% Get the correct figure dimensions (so complicated!)
h.screen = get(0,'screensize'); %How big is the screen?
h.figuresize = h.screen;
h.figuresize(3) = 0.75*h.screen(3); % Set figure width
h.figuresize(4) = 0.75*h.screen(4); % Set figure height.
h.fig = figure('Color',[1 1 1],...
                 'OuterPosition',h.figuresize,...
                 'Units','normalized',...
                 'Visible','off',...
                 'NumberTitle', 'off', ...
                 'Name','Slice viewer'); % Finally make the figure

movegui(h.fig,'center');
set(h.fig,'Visible','on'); %Make the figure visible after moving

% Make the main layout (using the GUI toolbox)

h.mainLayout = uix.VBoxFlex('Parent', h.fig, 'Spacing', 3);

% h.boxImage = uix.BoxPanel('Parent', h.mainLayout, 'Title', 'Slice');
% h.boxSelector = uix.BoxPanel('Parent', h.mainLayout, 'Title', 'Slice selector');

h.axImage = axes('Parent', uicontainer('Parent', h.mainLayout), ...
    'ActivePositionProperty', 'outerposition');

h.uiSlider = uicontrol('Parent', uicontainer('Parent', h.mainLayout), ...
    'Style', 'slider', ...
    'Tag', 'sliceSelector', ...
    'Units', 'normalized', ...
    'Callback', @onSliderChange, ...
    'Value', 1, ...
    'Max', size(ESI, 1), ...
    'Min', 1, ...
    'SliderStep', [1/(size(ESI, 1) - 1) 1/(size(ESI, 1) - 1)], ...
    'Position', [0.05 0.5 0.9 1]);

h.im = imagesc(h.axImage, squeeze(ESI(1, :, :)));
set(h.axImage, 'xtick', [], 'ytick', [], 'box', 'off', 'fontsize', 18);
colormap(h.axImage, 'gray');
axis(h.axImage, 'image');

set(h.mainLayout, 'Heights', [-15 -1]);

        % Add listeners
        try    % R2013b and older
           h.listeners.sliderVal = addlistener(h.uiSlider,...
                                     'ActionEvent',@onSliderChange);
           h.listeners.sliderVal = addlistener(h.uiSlider,...
                                     'ActionEvent',@onSliderChange);
        catch  % R2014a and newer
           h.listeners.sliderVal = addlistener(h.uiSlider,...
                                     'ContinuousValueChange',@onSliderChange);
           h.listeners.sliderVal = addlistener(h.uiSlider,...
                                     'ContinuousValueChange',@onSliderChange);
        end %Listeners

    function onSliderChange(source, ~)
%         h.axImage.Title = sprintf('Slice index %04d', get(source, 'Value'));
        index = round(get(source, 'Value'));
        title(h.axImage, sprintf('Slice index %04d', index));
        h.im.CData = squeeze(ESI(index, :, :));
    end

end


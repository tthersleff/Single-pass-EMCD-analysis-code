function [DP, ts] = getDiffImage(tagfilepath)
%GITDIFFIMAGE Imports a 4D STEM dataset into Matlab
%   gitDiffImage is a simple script written to import a 4D STEM dataset
%   into Matlab.  It uses timestamps that are recorded using code by Linus
%   Sch?nstr?m via Digital Micrograph hook-up scripts.  These timestamps
%   are written to a Gatan persistent tag structure that is saved on the
%   hard drive as a file with the *.gtg extention.
%   (c) 2018 Thomas Thersleff, Stockholm University

%% Import the tag structure
%   Get the tag file path
if ~exist('tagfilepath', 'var')
    [fName, fPath] = uigetfile({'*.gtg', 'Gatan Tag File'; ...
                                '*.*', 'All Files'}, ...
                                'Select the tag file');
    
    fString = strcat(fPath, fName);
else
    fString = tagfilepath;
end

fString = char(fString);

%   Get the tag structure and clean up
tags = dmread(fString);
clear fString fName fPath tagfile

%% Parse the tag file

% tags = tagfilepath; %Debug

TGfields = fieldnames(tags);
nTGfields = length(TGfields);

for idx = 1:nTGfields
    
    subTGfields = fieldnames(tags.(TGfields{idx}));
    subTGlength = length(subTGfields);
    if subTGlength < 1 
        T.(TGfields{idx}) = [];
    elseif subTGlength == 1 && strcmpi(subTGfields{1}, 'Value')
%         T{idx, 1} = tags.(TGfields{idx}).(subTGfields{1});
        T.(TGfields{idx}) = tags.(TGfields{idx}).(subTGfields{1});
    else
        temp = zeros(subTGlength, 1);
        for ind = 1:subTGlength
%             T{idx, ind} = tags.(TGfields{idx}).(subTGfields{ind}).Value;
            temp(ind) = tags.(TGfields{idx}).(subTGfields{ind}).Value;
        end
        T.(TGfields{idx}) = temp;
    end
    
    if iscell(T.(TGfields{idx}))
        T.(TGfields{idx}) = cell2mat(T.(TGfields{idx}));
    end
    
end

clear temp TGfields nTGfields ind idx


%   Find the first pixel index where the pixel beginning of the subsequent
%   pixel is less than the pixel end of the previous one.
extraPX = 1;
while mean(T.PixelBeginningTimestamps(extraPX:extraPX+9) - T.PixelEndTimestamps(1:10)) <= 0
    extraPX = extraPX + 1; 
end
extraPX = extraPX - 2;

if extraPX > 0
    T.PixelBeginningTimestamps(extraPX) = [];
elseif extraPX < 0
    T.PixelEndTimestamps(abs(extraPX)) = [];
end



tickres = 3330096; %ticks / s, computer dependent
firstpixel = T.PixelBeginningTimestamps(1);
ts.pxbegin = (T.PixelBeginningTimestamps - firstpixel) ./ tickres;
ts.pxend = (T.PixelEndTimestamps - firstpixel) ./ tickres;
ts.framegrab = (T.ImageRecorderTimestamps - firstpixel) ./ tickres;
ts.SIbegin = [0 (T.SIRestartTimestamps - firstpixel) ./ tickres]; %Number of times scanning started
ts.SIend = ([T.SIPauseTimestamps T.SIEndTimestamp] - firstpixel)./ tickres;
ts.SIpause = (T.SIPauseTimestamps - firstpixel) ./ tickres;
ts.SIresume = (T.SIRestartTimestamps - firstpixel) ./ tickres;
ts.startingPixels = ts.pxbegin(1); %Load with the first pixel of the entire ESI
nPix = double(prod(T.SurveyResolution));
nCols = double(T.SurveyResolution(2));
nRows = double(T.SurveyResolution(1));

numBadPx = 0;
if ~isempty(ts.SIpause)
    
    numBadPx = 0;
    for idx = 1:length(ts.SIpause)
        pxbegintspause = find(ts.pxbegin > ts.SIpause(idx), 1, 'first');
        pxendtspause = find(ts.pxend > ts.SIpause(idx), 1, 'first');
        remPX = logical(pxbegintspause - pxendtspause);
        if remPX
            ts.pxbegin(pxendtspause) = [];
            numBadPx = numBadPx + 1;
        end
        ts.startingPixels(idx + 1) = ts.pxbegin(pxendtspause);
    end
end

fprintf('Found and removed %d bad pixel timestamp(s)\n', extraPX + numBadPx);

% clear pxbegintspause pxendtspause didPauseAfterPXstarted;

%% Load the DP stack(s)

fprintf('ESI was paused %d time(s).  Expecting %d DP stack(s)\n', ...
    length(T.SIPauseTimestamps), length(ts.SIbegin));

for idx = 1:length(ts.SIbegin)
    [fName, fPath] = uigetfile({'*.dm4;*.dm3', 'Digital Micrograph file'; ...
                                '*.gtg', 'Gatan Tag File'; ...
                                '*.*', 'All Files'}, ...
                                sprintf('Select DP stack %d', idx));
    Atags{idx} = dmread(strcat(fPath, fName));

    qx = double(Atags{idx}.ImageList.Unnamed1.ImageData.Dimensions.Unnamed0.Value);
    qy = double(Atags{idx}.ImageList.Unnamed1.ImageData.Dimensions.Unnamed1.Value);
    nFrames(idx) = double(Atags{idx}.ImageList.Unnamed1.ImageData.Dimensions.Unnamed2.Value);
    
    Atemp = reshape(Atags{idx}.ImageList.Unnamed1.ImageData.Data.Value, ...
                qx, qy, nFrames(idx));
    
    % Ugly but functional way to find which frame was the last one acquired
    % in the frame grabber.  Useful if the view camera is stopped after
    % pausing the ESI acquisition, leaving empty frames but no time stamps.
    t = find(squeeze(sum(sum(Atemp, 2), 1)) ~= 0, 1, 'last');
    if isempty(t)
        finalframegrab(idx) = nFrames(idx);
    else
        finalframegrab(idx) = t;
    end
    
    if idx > 1
        A = cat(3, A, Atemp);
    else
        A = Atemp;
    end
end

clear Atemp fName fPath Atags;

%% Get the timestamp mask for the frame grabber

activeTS = false(size(ts.framegrab));

for idx = 1:length(ts.SIbegin)
    activePX(1:nFrames(idx), idx) = false(nFrames(idx), 1);
    firstframe(idx) = find(ts.framegrab > ts.startingPixels(idx), 1, 'first');
    lastframe(idx) = find(ts.framegrab < ts.SIend(idx), 1, 'last');
    activeTS(firstframe(idx) : lastframe(idx)) = true;
    if idx ~= 1
        firstframe(idx) = firstframe(idx) - finalframegrab(idx-1);
        lastframe(idx) = lastframe(idx) - finalframegrab(idx-1);
%         firstframe = firstframe - sum(nFrames(1: (idx - 1)));
%         lastframe = lastframe - sum(nFrames(1: (idx - 1)));
    end
    activePX(firstframe(idx) : lastframe(idx), idx) = true;
end

activePX = reshape(activePX, sum(nFrames), 1);

A(:, :, ~activePX) = [];
ts.activeframes = ts.framegrab(activeTS);
% for idx = 1:
% activePX = ts.framegrab >= 



% tags = dmread(tagfilepath);
% Atags = dmread(DPfile);
% A = reshape(Atags.ImageList.Unnamed1.ImageData.Data.Value, ...
%     Atags.ImageList.Unnamed1.ImageData.Dimensions.Unnamed0.Value, ...
%     Atags.ImageList.Unnamed1.ImageData.Dimensions.Unnamed1.Value, ...
%     Atags.ImageList.Unnamed1.ImageData.Dimensions.Unnamed2.Value);
% 
% fieldnames_pxbegin = fieldnames(tags.PixelBeginningTimestamps);
% fieldnames_pxend = fieldnames(tags.PixelEndTimestamps);
% fieldnames_framegrab = fieldnames(tags.ImageRecorderTimestamps);
% 
% % Ensure that the fields have the same length
% if length(fieldnames_pxbegin) ~= length(fieldnames_pxend)
%     fieldnames_pxbegin(1) = []; % Remove the first field
% end
% 
% for idx = 1:length(fieldnames_pxbegin)
%     pxbegin(idx) = tags.PixelBeginningTimestamps.(fieldnames_pxbegin{idx}).Value;
%     pxend(idx) = tags.PixelEndTimestamps.(fieldnames_pxend{idx}).Value;
% end
% 
% for idx = 1:length(fieldnames_framegrab)
%     framegrab(idx) = tags.ImageRecorderTimestamps.(fieldnames_framegrab{idx}).Value; 
% end
% 
% nCols = dims(2);
% nRows = dims(1);
% nPix = prod(dims);
% [qx, qy, ~] = size(A);


ts.finalrowindex = 1 : nCols : nPix+1;
ts.timestampNewRow = ts.pxbegin(ts.finalrowindex(1:end-1));

% lastframe = find(framegrab > pxend(end), 1, 'first');

for idx = 1:length(ts.timestampNewRow)
    ts.idxNewRow(idx) = find(ts.activeframes >= ts.timestampNewRow(idx), 1, 'first'); 
end

% %   If the frame grabber finished before the pixel scanner, then the last
% %   frame needs to be the last frame of the framegrabber.  
% if isempty(lastframe)
%     lastframe = length(framegrab) - 1;
% end

ts.idxNewRow = ts.idxNewRow + 1;
% ts.idxNewRow(end+1) = lastframe + 1;
ts.idxNewRow(end+1) = length(ts.activeframes) + 1;

DP = zeros(qx, qy, nPix);
ts.maphistogram = zeros(nCols, nRows);

for idx = 1:nRows
    temp = A(:, :, ts.idxNewRow(idx):(ts.idxNewRow(idx+1)-1));
    temp2 = resampleESI(temp, [qx qy nCols], 'nearest');
    DP(:, :, ts.finalrowindex(idx):ts.finalrowindex(idx+1)-1) = temp2;
    ts.maphistogram(:, idx) = ...
        histcounts(ts.activeframes(ts.idxNewRow(idx) : (ts.idxNewRow(idx+1)-1)), ...
        [ts.pxbegin(ts.finalrowindex(idx) : ts.finalrowindex(idx+1)-1); ts.pxend(ts.finalrowindex(idx+1)-1)]);
end

DP = reshape(DP, qx, qy, nCols, nRows);

% %% Plot figure
% 
% figure('Color', 'w', 'Position', [271,199,1071,731]);
% plot(1:nPix, pxbegin, 'k.', 'DisplayName', 'Pixel Begin');
% hold on;
% plot(1:nPix, pxend, 'r.', 'DisplayName', 'Pixel End');
% plot(1:length(framegrab), framegrab, 'bo', 'DisplayName', 'Image Recorder');
% title('Time stamp comparison');
% xlabel('Linear index');
% ylabel('Timestamp');
% set(gca, 'fontsize', 18);
% grid on;
% 
% figure('Color', 'w', 'Position', [271,199,1071,731]);
% plot(1:nPix, pxend - pxbegin, 'k.', 'DisplayName', 'Pixel Dwell Time');
% title('Pixel Dwell Time');
% xlabel('Linear index');
% ylabel('Dwell time (timestamp units)');
% set(gca, 'fontsize', 18);
% grid on;

end


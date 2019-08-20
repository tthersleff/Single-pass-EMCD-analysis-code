function [v] = makeMovie(filmStruct, movieName, frameRate)
%MAKEMOVIE Simple code to make a movie from a film struct

v = VideoWriter(movieName, 'MPEG-4');

if ~exist('frameRate', 'var')
    v.FrameRate = 20;
else
    v.FrameRate = frameRate;
end

v.Quality = 20;

open(v);
writeVideo(v, filmStruct);
close(v);

end


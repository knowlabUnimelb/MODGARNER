% SHOWTEXT Display text
%       To be used with Psychophysics Toolbox .mex file, SCREEN
%       SHOWTEXT (screenparms, size, stimulus, onoff, vertdown, leftright)
%       displays 'stimulus' in a window defined in 'screenparms' where
%       'size' = stimulus size in pixels
%
%       'onoff' specifies color (0 = white, 1 = black)
%       'vertdown' specifies vertical offset from center in pixels
%       'leftright' specifies horizontal offset from center in pixels
%
%       See also prepexp, screen

function showtext (screenparms, size, stimulus, varargin)

optargs = {1, 0, 0, screenparms.sansSerifFont};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[color, vertdown, leftright, font] = optargs{:}; % Place optional args in memorable variable names

%%
size = round(size);
if color == 1
    c = screenparms.black;
elseif color == 0
    c = screenparms.white;
else
    c = color;
end

fromtop = (screenparms.rect(RectBottom) - screenparms.rect(RectTop))/2 + vertdown;
fromleft = (screenparms.rect(RectRight) - screenparms.rect(RectLeft))/2 - (length(stimulus)/2 * size/2.4) + leftright;

Screen('TextSize', screenparms.window, size);
Screen('TextFont', screenparms.window, font);
Screen('DrawText', screenparms.window, stimulus, fromleft, fromtop, c);

function [h, probs, binsc] = plotHist(x, varargin)

optargs = {.02, 0, 1, '-r', 2, 1};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin
[eps, minx, maxx, linestr, lw, wp] = optargs{:}; % Place optional args in memorable variable names

%%
binse = minx:eps:maxx;     % Bin edges
binsc = minx+eps/2:eps:maxx-eps/2; % Bin centres

count = histc(x,binse); % Count samples in between bins
count = count(1:end-1);                           % Drop last bin (see help histc)
probs = count/sum(count)/eps;                     % Convert bin frequencies to probabilities

if wp == 1
    h = plot(binsc, probs, linestr, 'LineWidth', lw);             % Plot probabilities at bin centres
else
    h = bar(binsc, probs, 'hist');
end

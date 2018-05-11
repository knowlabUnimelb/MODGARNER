function [data, cols, mds] = getPosteriorPredictionTrialData(dataloc, datafn, mdsloc, mdsfn, reversexy)
% Data loc is string containing data folder location
% Data fn is the name of the data file to load

data = dlmread(fullfile(dataloc, datafn)); % Load data
data(1,:) = []; % First row is a count of the number of columns (for Nosofsky FORTRAN, delete here)
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc', 'igp'}; % Data columns

%% Organize data columns
item = data(:, strcmp(cols, 'item'));     % Item
rt   = data(:, strcmp(cols, 'rts'))/1000; % RT

nD   = size(data, 1);% Number of data points
cat  = ones(nD, 1);  % Category for each item
cat(item >= 7) = 2;  % Category for B items equals 2
opp = 2 - (cat - 1); % Opposite category for each item

% Set up output
cols = [cols, 'cat', 'opp'];
data(:,strcmp(cols, 'rts')) = rt;
data(:,strcmp(cols, 'cat')) = cat;
data(:,strcmp(cols, 'opp')) = opp;

%% Load MDS
% Note that in Garner Within - Saturation condition, the stimulus
% coordinates need to be reversed for X and Y as follows:
% stimulusCoordinates = stimulusCoordinates(:, [2 1]);
stimulusCoordinates = load(fullfile(mdsloc, mdsfn));
if reversexy
    mds = stimulusCoordinates(:, [2 1]);
else
    mds = stimulusCoordinates;
end
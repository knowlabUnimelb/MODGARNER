% Load data from Garner Within (Integral) Experiment
% Condition: B = Brightness, S = Saturation
%

function [data, cols] = readGarnerExperimentData(folder, dataPrefix, subStr, subject_number)
% Subject number refers to an entry in subs

%% Subject data file details
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc'}; % Data columns
% cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc', 'igp'}; % Data columns
% The igp column flags the first trial in every task but it is not a column
% in the datafile so don't add it yet

%% Load data
sub = subStr{subject_number};

fn = sprintf(dataPrefix, sub); % File name sequence data, tagged with outliers
data = dlmread(fullfile(folder, fn));                      % Load data
data(1,:) = [];                                            % Remove nLines number at top of file
% The first line of the datafiles is a number of lines (used in Nosofsky's FORTRAN code but not needed here)

ncols = size(data, 2);
data(:, ncols + 1) = zeros(size(data, 1), 1);              % Preallocate
data(data(:,strcmp(cols, 'trl')) == 1, ncols+1) = 1;       % Set first trial to 1 to ignore
cols = [cols, 'igp'];                                      % Add ignore column                                 
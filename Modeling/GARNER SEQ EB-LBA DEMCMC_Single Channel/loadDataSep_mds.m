% Load separable data
function [model, data, n, subnums] = loadDataSep

model = @seqEBLBA;

cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc'};
% subs = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'};
% subnums = [101:104 201:204];

subs = {'B1', 'B2', 'B3', 'B4'};
subnums = [101:104];

% subs = {'S1', 'S2', 'S3', 'S4'};
% subnums = [201:204];

data = struct('data', [], 'cols', []);
for sidx = 1:numel(subs)
    % Load Data information
    subject = subs{sidx};
    
    fn = sprintf('Garner_SEP_sub%s_ordered_tagged_msec.dat', subject); % Data filename
    sdata = dlmread(fullfile(pwd, 'Data', fn));
    sdata(1,:) = []; % Remove nLines number at top of file
    sdata(:, strcmp(cols, 'rts'))   = sdata(:, strcmp(cols, 'rts'))./1000; % Convert RT to sec
    
    igp  = sdata(:, strcmp(cols, 'trl')) == 1; % Equals 1 if it's the start of a new block
    cat  = ones(size(sdata, 1), 1);  % Category for each item
    cat(sdata(:,strcmp(cols, 'item')) >= 7) = 2;
    opp = 2 - (cat - 1);             % Opposite category for each item
    
    sdata = [sdata, igp, cat, opp];
    
    data.data = [data.data; sdata];
    
    %% MDS information
%     if strcmp(subject(1), 'B')
        data.stimloc{sidx} = load(fullfile(pwd, 'MDS', 'mds_Separable.dat'));
%     else % If saturation, flip xy
        data.stimloc{sidx} = load(fullfile(pwd, 'MDS', 'mds_Separable.dat'));
%         data.stimloc{sidx} = data.stimloc{sidx}(:, [2 1]);
%     end
end
cols = [cols, 'igp', 'cat', 'opp'];
data.cols = cols;
data.model = model;

n.subjects = numel(subs);
n.items = size(data.stimloc{1},1);

%% Mean center each subject at the grand mean
rtmeans = aggregate(data.data(data.data(:,strcmp(data.cols, 'exc')) ~= 1, :),...
    mstrfind(data.cols, {'sub'}), mstrfind(data.cols, {'rts'}), @mean, 1);
rtcnts  = aggregate(data.data,...
    mstrfind(data.cols, {'sub'}), mstrfind(data.cols, {'rts'}), @count, 1);

submeans = [];
for i = 1:numel(rtmeans)
    submeans = [submeans; repmat(rtmeans(i), rtcnts(i), 1)];
end

overallmean = mean(data.data(data.data(:,strcmp(data.cols, 'exc')) ~= 1, strcmp(data.cols, 'rts')));

data.data(:, strcmp(data.cols, 'rts')) = data.data(:, strcmp(data.cols, 'rts')) - submeans + overallmean;
    
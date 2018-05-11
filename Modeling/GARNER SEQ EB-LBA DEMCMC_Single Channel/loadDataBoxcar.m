% Load separable data
function [model, data, n, subnums] = loadDataBoxcar

model = @seqEBLBA_single;

cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc'};
subs = {'L1', 'L2', 'L3', 'L4', 'S1', 'S2', 'S3', 'S4'};
subnums = [101:104 205 202:204];

% subs = {'L1', 'L2', 'L3', 'L4'};
% subnums = [101:104];

% subs = {'S1', 'S2', 'S3', 'S4'};
% subnums = [205 202:204];


data = struct('data', [], 'cols', []);
for sidx = 1:numel(subs)
    % Load Data information
    subject = subs{sidx};
    
    fn = sprintf('Garner_BOXCAR_sub%s_ordered_tagged_msec.dat', subject); % Data filename
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
%     if strcmp(subject(1), 'L')
        data.stimloc{sidx} = load(fullfile(pwd, 'MDS', 'mds_Separable.dat'));
%     else % If saturation, flip xy
%         data.stimloc{sidx} = load(fullfile(pwd, 'MDS', 'mds_Separable.dat'));
%         data.stimloc{sidx} = data.stimloc{sidx}(:, [2 1]);
%     end
end
cols = [cols, 'igp', 'cat', 'opp'];
data.cols = cols;
data.model = model;

n.subjects = numel(subs);
n.items = size(data.stimloc{1},1);

%% Fix subnum issue so that rtmeans are ordered properly for mean centering
if strcmp(subs{1}, 'S1')
data.data(data.data(:,strcmp(cols, 'sub')) == 205, strcmp(cols, 'sub')) = 201;
subnums = [201:204];
elseif strcmp(subs{5}, 'S1')
data.data(data.data(:,strcmp(cols, 'sub')) == 205, strcmp(cols, 'sub')) = 201;
subnums = [101:104  201:204];
end


%% Mean center each subject at the grand mean
rtmeans = aggregate(data.data(data.data(:,strcmp(data.cols, 'exc')) ~= 1, :),...
    mstrfind(data.cols, {'sub'}), mstrfind(data.cols, {'rts'}), @nanmean, 1);
rtcnts  = aggregate(data.data,...
    mstrfind(data.cols, {'sub'}), mstrfind(data.cols, {'rts'}), @count, 1);

submeans = [];
for i = 1:numel(rtmeans)
    submeans = [submeans; repmat(rtmeans(i), rtcnts(i), 1)];
end

overallmean = nanmean(data.data(data.data(:,strcmp(data.cols, 'exc')) ~= 1, strcmp(data.cols, 'rts')));

data.data(:, strcmp(data.cols, 'rts')) = data.data(:, strcmp(data.cols, 'rts')) - submeans + overallmean;

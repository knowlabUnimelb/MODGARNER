function [model, data, n, subnums] = loadData

model = @seqEBLBA;
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc'};
condition = 'all'; % brightness, saturation, all
mdsVer = 'Constrained'; % {'Full', 'Constrained'}; % Note: Use Constrained

switch condition
    case 'all'
        subs = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'};
        subnums = [101 103 105 107 102 104 106 108];
    case 'brightness'
        subs = {'B1', 'B2', 'B3', 'B4'};
        subnums = 101:2:107;
    case 'saturation'        
        subs = {'S1', 'S2', 'S3', 'S4'};
        subnums = 102:2:108;
end

%% Load data
data = struct('data', [], 'cols', []);
for sidx = 1:numel(subs)
    % Load Data information
    subject = subs{sidx};
    
    fn = sprintf('Garner_Exp1_sub%s_ordered_tagged.dat', subject); % Data filename
    sdata = dlmread(fullfile(pwd, 'Data', fn));
    sdata(1,:) = []; % Remove nLines number at top of file (useful only for Nosofsky FORTRAN)
    sdata(:, strcmp(cols, 'rts'))   = sdata(:, strcmp(cols, 'rts'))/1000; % Convert RT to sec
    
    igp  = sdata(:, strcmp(cols, 'trl')) == 1; % Equals 1 if it's the start of a new block
    cat  = ones(size(sdata, 1), 1);  % Category for each item
    cat(sdata(:,strcmp(cols, 'item')) >= 7) = 2;
    opp = 2 - (cat - 1);             % Opposite category for each item
    
    sdata = [sdata, igp, cat, opp];
    
    data.data = [data.data; sdata];
    
    %% MDS information
    if strcmp(subject(1), 'B')
        data.stimloc{sidx} = load(fullfile(pwd, 'MDS', sprintf('mds_Brightness - %s.dat', mdsVer)));
    else % If saturation, flip xy 
        data.stimloc{sidx} = load(fullfile(pwd, 'MDS', sprintf('mds_Saturation - %s.dat', mdsVer)));
        data.stimloc{sidx} = data.stimloc{sidx}(:, [2 1]);
    end
end
cols = [cols, 'igp', 'cat', 'opp'];
data.cols = cols;
data.model = model;

n.subjects = numel(subs);
n.items = size(data.stimloc{1},1);

if numel(subs) == 8
    data.data(ismember(data.data(:,strcmp(cols, 'sub')), [102 104 106 108]), strcmp(cols, 'sub')) =...
        data.data(ismember(data.data(:,strcmp(cols, 'sub')), [102 104 106 108]), strcmp(cols, 'sub')) + 100;
    subnums = [101 103 105 107 202 204 206 208];
end

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
    
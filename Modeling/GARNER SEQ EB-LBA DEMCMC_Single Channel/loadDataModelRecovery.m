function [model, data, n, subnums] = loadDataModelRecovery

model = @seqEBLBA;
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc'};
condition = 'simulation'; % brightness, saturation, all
mdsVer = 'Constrained'; % {'Full', 'Constrained'}; % Note: Use Constrained

switch condition
    case 'simulation'
        subs = {'S01', 'S02', 'S03', 'S04', 'S05'};
        subnums = [201:205];
    
end

%% Load data
data = struct('data', [], 'cols', []);
for sidx = 1:numel(subs)
    % Load Data information
    subject = subs{sidx};
    
    fn = sprintf('Garner_SIM2_sub%s_ordered_tagged.dat', subject); % Data filename
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
            data.stimloc{sidx} = load(fullfile(pwd, 'MDS', 'mds_Separable.dat'));
%     if strcmp(subject(1), 'B')
%         data.stimloc{sidx} = load(fullfile(pwd, 'MDS', sprintf('mds_Brightness - %s.dat', mdsVer)));
%     else % If saturation, flip xy 
%         data.stimloc{sidx} = load(fullfile(pwd, 'MDS', sprintf('mds_Saturation - %s.dat', mdsVer)));
%         data.stimloc{sidx} = data.stimloc{sidx}(:, [2 1]);
%     end
end
cols = [cols, 'igp', 'cat', 'opp'];
data.cols = cols;
data.model = model;

n.subjects = numel(subs);
n.items = size(data.stimloc{1},1);
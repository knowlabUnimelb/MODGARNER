function [means, accuracy] = getSubjectPredictionAverages(data, cols, simrt, simresp)

nPosteriorPredictionSamples = size(simrt, 2);

% Build column labels for simrt and simresp
labelRT  = @(x)(sprintf('simRT%03d', x));
labelACC = @(x)(sprintf('simACC%03d', x));

simRTlabels = cellfun(labelRT, mat2cell(1:nPosteriorPredictionSamples, 1, ones(1, nPosteriorPredictionSamples)), 'UniformOutput', false);
simACClabels = cellfun(labelACC, mat2cell(1:nPosteriorPredictionSamples, 1, ones(1, nPosteriorPredictionSamples)), 'UniformOutput', false);

cols = [cols, simRTlabels, simACClabels];
data = [data, simrt, simresp];

%% Separate data into cells for different tasks
ntasks = numel(unique(data(:, mstrfind(cols, {'task'}))));
sessTask = unique(data(:, mstrfind(cols, {'sess', 'task'})),'rows'); % Indexes for each session and task

% Note: original code used the following columns in taskData
%   {'trl', 'item', 'acc', 'rts', 'exc', 'igp'}
taskData = cell(1,ntasks); % Preallocate
cols = [cols, 'prev']; % Update cols for newTaskData cell matrix below
for i = 1:ntasks
    % Allocate each task to a cell
    taskData(i) = mat2cell(data(data(:,strcmp(cols, 'task')) == i, :),...
        size(data(data(:,strcmp(cols, 'task')) == i, :), 1),...
        size(data(data(:,strcmp(cols, 'task')) == i, :), 2));
    
    % Recode items into previous and current
    [newTaskData{i}, cols] = codePreviousTrial(taskData{i}, cols);
    
    % Recode control2 --> control1, correlated2 --> correlated 1, and
    % align irrelevant dimension change in the filtering condition
    newTaskData{i} = recodeTaskItems(newTaskData{i}, cols);
end

% Combine control2 with control1 and correlated2 with correlated1
finalTaskData{1} = [newTaskData{1}; newTaskData{2}];
finalTaskData{2} = [newTaskData{3}; newTaskData{4}];
finalTaskData{3} = [newTaskData{5}];
nTaskConditions  = size(finalTaskData,2);

for i = 1:nTaskConditions
    %     finalTaskData = newTaskData;
    %     for i = 1:ntasks
    % Average errors after removing excluded trials
    accuracy{i} =...
        aggregate(finalTaskData{i}(all(finalTaskData{i}(:,mstrfind(cols, {'exc', 'igp'})) == 0, 2), :),...
        mstrfind(cols, {'prev', 'item'}),...
        find(strncmp(cols, 'simACC', 6)),...
        @(x)(nanmean(1-x)));
        
    % Average rts after removing excluded trials and errors
    means{i} =...
        aggregate(finalTaskData{i}(all([finalTaskData{i}(:,mstrfind(cols, {'exc', 'igp'})) == 0, finalTaskData{i}(:,strcmp(cols, 'acc')) == 1], 2), :),...
        mstrfind(cols, {'prev', 'item'}),...
        find(strncmp(cols, 'simRT', 5)),...
        @nanmean);
end

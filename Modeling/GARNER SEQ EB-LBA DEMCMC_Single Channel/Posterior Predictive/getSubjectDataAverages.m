function [means, accuracy, subcheck] = getSubjectDataAverages(datafolder, dataPrefix, subStr, subIdxs)
% clear all
% clc
% close all

% mainfolder = 'C:\Users\littled\Dropbox\Work\2014 Garner\GarnerDEMCMCsep';
% datafolder = fullfile(mainfolder, 'Data');
% subIdxs = 5:8; % 1:4 - Brightness, 5:8 - Saturation
% outputFileName = 'GarnerWithin_Sat';

%% %%%%%%%%%%%%%%%%%%%%%%%%% ORGANIZE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle through each subject (if more than one indicated)
for sidx = 1:numel(subIdxs);
    % Read raw data
    [data, cols] = readGarnerExperimentData(datafolder, dataPrefix, subStr, subIdxs(sidx));
    
    subcheck(sidx) = data(1, strcmp(cols, 'sub'));
    
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
        accuracy{i}(:,:, sidx) =...
            aggregate(finalTaskData{i}(all(finalTaskData{i}(:,mstrfind(cols, {'exc', 'igp'})) == 0, 2), :),...
            mstrfind(cols, {'prev', 'item'}),...
            strcmp(cols, 'acc'),...
            @(x)(mean(1-x)));
        
        % Average rts after removing excluded trials and errors
        means{i}(:,:,sidx) =...
            aggregate(finalTaskData{i}(all([finalTaskData{i}(:,mstrfind(cols, {'exc', 'igp'})) == 0, finalTaskData{i}(:,strcmp(cols, 'acc')) == 1], 2), :),...
            mstrfind(cols, {'prev', 'item'}),...
            strcmp(cols, 'rts'),...
            @mean);
    end
end
% Note: Combining the recoded data and then averaging produces slightly
% different numbers compared to averaging, recoding, then averaging across
% tasks due to differences in the trial numbers

% %% Average across subjects
% favgSub = @(x)(mean(x, 3));
% fstdSub = @(x)(std(x, [], 3));
% 
% ma = cellfun(favgSub, accuracy, 'UniformOutput', false); % Accuracy 
% sa = cellfun(fstdSub, accuracy, 'UniformOutput', false); % Accuracy std
% ca = mat2cell(numel(subIdxs) * ones(1,nTaskConditions), 1, ones(1, nTaskConditions)); % Count
% mm = cellfun(favgSub, means, 'UniformOutput', false); % RT 
% sm = cellfun(fstdSub, means, 'UniformOutput', false); % RT std 
% cm = ca; 

% 
% %% Save output
% save(outputFileName)
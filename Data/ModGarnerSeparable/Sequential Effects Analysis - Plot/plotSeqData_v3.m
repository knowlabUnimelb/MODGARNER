% Plot Garner Data
% Average over different conditions
clear all
clc
close all

plotQuantile = false;
plotmeans = true;
plotaccuracy = true;

%% Initial Sequp
subs =  {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'};
itemSets = {1:2:12; 2:2:12; [1 3 5 8 10 12]; [2 4 6 7 9 11]; 1:12}; % Sets of items for each condition
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc', 'igp'}; % Data columns

% subIdxs = 1:4;
subIdxs = 5:8;        
dataloc = 'C:\Users\littled\Dropbox\work\2017 Garner\MODIFIED GARNER SEPARABLE\data';

%% Get Grand mean
alldata =[];
for sidx = subIdxs
    sub = subs{sidx};
    fn = sprintf('Garner_SEP_sub%s_ordered_tagged_msec.dat', sub); % File name sequence data, tagged with outliers
    data = dlmread(fullfile(dataloc, fn));                 % Load data
    data(1,:) = [];                                            % Remove nLines number at top of file
    
    alldata = [alldata; data];
end
grandmean = nanmean(alldata(alldata(:,end)~=1,end-1));

%% Cycle through each subject
tic
scnt = 1;
for sidx = subIdxs
    sub = subs{sidx};
    
    %% Load data
    fn = sprintf('Garner_SEP_sub%s_ordered_tagged_msec.dat', sub); % File name sequence data, tagged with outliers
    data = dlmread(fullfile(dataloc, fn));                 % Load data
    data(1,:) = [];                                            % Remove nLines number at top of file
   
    subjectmean(scnt) = nanmean(data(data(:,strcmp(cols, 'exc')) ~= 1, strcmp(cols, 'rts')));
%     data(:, strcmp(cols, 'rts')) = data(:, strcmp(cols, 'rts')) - nanmean(data(:, strcmp(cols, 'rts'))) + (grandmean * 1000); % grandmean = .7525;
    
    %% Organize data columns
    exc  = data(:,strcmp(cols, 'exc'));       % Outlying RT tag
    igp  = data(:, strcmp(cols, 'trl')) == 1; % Equals 1 if it's the start of a new block
    task = data(:, strcmp(cols, 'task'));     % Task
    item = data(:, strcmp(cols, 'item'));     % Item
    acc  = data(:, strcmp(cols, 'acc'));      % Accuracy
    rt   = data(:, strcmp(cols, 'rts'))/1000; % RT
    
    nD   = size(data, 1);% Number of data points
    cat  = ones(nD, 1);  % Category for each item
    cat(item >= 7) = 2;  % Category for B items equals 2
    opp = 2 - (cat - 1); % Opposite category for each item
    
    %% Separate data into cells for different tasks
    sessTask = unique(data(:, mstrfind(cols, {'sess', 'task'})),'rows'); % Indexes for each session and task
    taskData = cell(1,5); % Preallocated
    
    % Cycle through each session and task
    for i = 1:size(sessTask, 1)
        if sessTask(i,1) == 1; % Preallocate empty matrix if the session equals 1
            taskData{sessTask(i,2)} = [];
        end
        
        % Extract data by session & task. Also record whether item should
        % be ignored because it is the first item of a block
        st = [data(data(:, strcmp(cols,'sess')) == sessTask(i,1) & data(:, strcmp(cols, 'task')) == sessTask(i,2), :),...
            igp(data(:, strcmp(cols,'sess')) == sessTask(i,1) & data(:, strcmp(cols, 'task')) == sessTask(i,2))];
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        st(1,mstrfind(cols, {'acc', 'rts'})) = nan(1,2); % Set first trial from each task and block to nan to ignore
        % Why? I guess because I can't plot the data by the previous item
        % if there is no previous item
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Add task/block trials to current task matrix
        taskData{sessTask(i,2)} = [taskData{sessTask(i,2)}; st(:, mstrfind(cols, {'trl', 'item', 'acc', 'rts', 'exc', 'igp'}))];
    end
    
    %% Recode B items and averaged
    % Cycle through tasks
    for i = 1:5
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        taskData{i}(taskData{i}(:,3) == 2, 3) = 0; % Recode B items from 2 to 0
        % Why? There's no category information in taskdata
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         taskData{i}(isnan(taskData{i}(:,4)),:) = []; % Remove nans
        %         taskData{i}(taskData{i}(:,4) == 0,:) = []; % Remove errors
        
        %% Flip all of the items so that I'm only plotting same category (on the left) and different category (on the right) rather than A or B
        prev = taskData{i}(1:end-1,2);        % Items on previous trial
        currData = taskData{i}(2:end, 2:end); % Items on current trial
        items = unique(prev);                 % Item identifiers
        if i < 5 % Task is control or correlated
            rev   = flipud(items); % Flip the items so that they match with their corresponding other category item
        else
            rev = [11 12 9 10 7 8 5 6 3 4 1 2]'; % This is more complex for the filtration condition since only one axes need to be flipped
        end
        
       
        % Recode items so that you have three current item types and all
        % previous items are reversed
        for cdidx = 1:size(currData,1)
            if currData(cdidx,1) >= 7
                prev(cdidx,1) = rev(prev(cdidx,1) == items);
                currData(cdidx,1) = rev(currData(cdidx,1) == items);
            end
        end
        
        prevItems = unique(prev);          % Unique identifiers for previous items
        iitems    = unique(currData(:,1)); % Unique identifiers for flipped items
        qcnt = 1;
        for pi = 1:numel(prevItems)
            for pj = 1:numel(iitems)
                % Compute mean 
                accuracy{i}(qcnt, :, scnt) = [prevItems(pi), iitems(pj),...
                    numel(currData(prev == prevItems(pi) & currData(:,1) == iitems(pj) &...
                    currData(:,end-1) == 0 & currData(:,end) == 0 & ~isnan(currData(:,2)) & currData(:,2) == 0, 3))./...
                    numel(currData(prev == prevItems(pi) & currData(:,1) == iitems(pj) &...
                    currData(:,end-1) == 0 & currData(:,end) == 0 & ~isnan(currData(:,2)), 3))];
                means{i}(qcnt, :, scnt) = [prevItems(pi), iitems(pj),...
                    mean(currData(prev == prevItems(pi) & currData(:,1) == iitems(pj) &...
                    currData(:,end-1) == 0 & currData(:,end) == 0 & ~isnan(currData(:,2)) & currData(:,2) == 1, 3))];
                quantiles{i}(qcnt, :, scnt) = [prevItems(pi), iitems(pj),...
                    prctile(currData(prev == prevItems(pi) & currData(:,1) == iitems(pj) &...
                    currData(:,end-1) == 0 & currData(:,end) == 0 & ~isnan(currData(:,2)) & currData(:,2) == 1, 3), [10:20:90])];
                qcnt = qcnt + 1;
            end
        end
    end
    
    toc
    scnt = scnt + 1;
end

%% Average across tasks
for i = 1:(scnt-1)
    % Accuracy
    newaccuracy{1}(:,:,i) = [accuracy{1}(:,1:2,i), mean([accuracy{1}(:,3,i), accuracy{2}(:,3,i)], 2)];
    newaccuracy{2}(:,:,i) = [accuracy{3}(:,1:2,i), mean([accuracy{3}(:,3,i), accuracy{4}(:,3,i)], 2)];
    
    temp = accuracy{5}(:,:,i);
    
    remap1 = [2 1; 4 3; 6 5];
    remap2 = [1 2; 2 1; 3 4; 4 3; 5 6; 6 5; 7 8; 8 7; 9 10; 10 9; 11 12; 12 11];
    remap3 = [1 1; 2 1; 3 3; 4 3; 5 5; 6 5; 7 7; 8 7; 9 9; 10 9; 11 11; 12 11];
    
    % Remap boundary items for filtration condition maintaining irrelevant
    % dimension
    for j = 1:size(remap2,1)
       temp( ismember(accuracy{5}(:,2,i), [6]) & ismember(accuracy{5}(:,1,i), remap2(j,1)), 1) = remap2(j,2); 
    end
    % Collapse irrelevant dimension for nonboundary items
    for j = 1:size(remap3,1)
       temp( ismember(accuracy{5}(:,2,i), [1:4]) & ismember(accuracy{5}(:,1,i), remap3(j,1)), 1) = remap3(j,2); 
    end
    
    for j = 1:size(remap1,1)
        temp(ismember(accuracy{5}(:,2,i), remap1(j,1)),2) = remap1(j,2);
    end
    temp = aggregate(temp, [1 2], 3);
    newaccuracy{3}(:,:,i) = temp;
    
    % Mean RTs
    newmeans{1}(:,:,i) = [means{1}(:,1:2,i), mean([means{1}(:,3,i), means{2}(:,3,i)], 2)];
    newmeans{2}(:,:,i) = [means{3}(:,1:2,i), mean([means{3}(:,3,i), means{4}(:,3,i)], 2)];
    
    temp = means{5}(:,:,i);
    
%     remap1 = [2 1; 4 3; 6 5];
%     remap2 = [1 2; 2 1; 3 4; 4 3; 5 6; 6 5; 7 8; 8 7; 9 10; 10 9; 11 12; 12 11];
    
    for j = 1:size(remap2,1)
       temp(ismember(means{5}(:,2,i), [6]) & ismember(means{5}(:,1,i), remap2(j,1)), 1) = remap2(j,2); 
    end
    for j = 1:size(remap3,1)
       temp( ismember(means{5}(:,2,i), [1:4]) & ismember(means{5}(:,1,i), remap3(j,1)), 1) = remap3(j,2); 
    end
    for j = 1:size(remap1,1)
        temp(ismember(means{5}(:,2,i), remap1(j,1)),2) = remap1(j,2);
    end
    temp = aggregate(temp, [1 2], 3);
    newmeans{3}(:,:,i) = temp;
    
    % Quantiles
    newquantiles{1}(:,:,i) = [quantiles{1}(:,1:2,i), [quantiles{1}(:,3:7,i) + quantiles{2}(:,3:7,i)]./2];
    newquantiles{2}(:,:,i) = [quantiles{3}(:,1:2,i), [quantiles{3}(:,3:7,i) + quantiles{4}(:,3:7,i)]./2];
    
    temp = quantiles{5}(:,:,i);
    
%     remap1 = [2 1; 4 3; 6 5];
%     remap2 = [1 2; 2 1; 3 4; 4 3; 5 6; 6 5; 7 8; 8 7; 9 10; 10 9; 11 12; 12 11];
    
    for j = 1:size(remap2,1)
       temp(ismember(quantiles{5}(:,2,i), [2 4 6]) & ismember(quantiles{5}(:,1,i), remap2(j,1)), 1) = remap2(j,2); 
    end
    for j = 1:size(remap1,1)
        temp(ismember(quantiles{5}(:,2,i), remap1(j,1)),2) = remap1(j,2);
    end
    temp = aggregate(temp, [1 2], 3:7);
    newquantiles{3}(:,:,i) = temp;
end

%% Average across subjects
for i = 1:3
    ma{i} = nanmean(newaccuracy{i}, 3);
    sa{i} = nanstd(newaccuracy{i}, [], 3);
    ca{i} = scnt - 1;
    
    mm{i} = nanmean(newmeans{i}, 3);
%     sm{i} = nanstd(newmeans{i}, [], 3);
    
    withinmeans{i} = newmeans{i};
    for j = 1:4
        withinmeans{i}(:,3,j) = newmeans{i}(:,3,j) - subjectmean(j) + grandmean;
    end
    sm{i} = nanstd(withinmeans{i}, [], 3);
    cm{i} = scnt - 1;
    
    mq{i} = nanmean(newquantiles{i}, 3);
    sq{i} = nanstd(newquantiles{i}, [], 3);
    cq{i} = scnt - 1;
end

%% Quantile plots
if plotQuantile
    % Control task
    figure('WindowStyle', 'docked')
    tasks = [1 1 1];
    items = [1 3 5];
    itemtitles = {'Far', 'Mid', 'Near'};
    for i = 1:3
        subplot(3,3,i)
        garnerQPplotv2(mq, tasks(i), items(i), '-ok')
%         hold on
%         garnerQPplotv2(pq, tasks(i), items(i), '-xk')
        set(gca,'YLim', [300 1500])
        set(gca,'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'})
        title(sprintf('%s', itemtitles{mod(i-1, 3)+1}), 'FontSize', 14)
    end
    
    % Correlated task
    tasks = [2 2 2];
    items = [1 3 5];
    itemtitles = {'Far', 'Mid', 'Near'};
    for i = 1:3
        subplot(3,3,i+3)
        garnerQPplotv2(mq, tasks(i), items(i), '-ok')
%         hold on
%         garnerQPplotv2(pq, tasks(i), items(i), '-xk')
        set(gca,'YLim', [300 1500])
        set(gca,'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'})
        title(sprintf('%s', itemtitles{mod(i-1, 3)+1}), 'FontSize', 14)
    end
    
    % Filtering task
    tasks = [3 3 3];
    items = [1 3 5];
    itemtitles = {'Far', 'Mid', 'Near'};
    for i = 1:3
        subplot(3,3,i+6)
        garnerQPplotv2(mq, tasks(i), items(i), {'-ok', '-or'})
%         hold on
%         garnerQPplotv2(pq, tasks(i), items(i), {'-xk', '-xr'})
        set(gca,'YLim', [300 1700])
        set(gca,'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'})
        title(sprintf('%s', itemtitles{mod(i-1, 3)+1}), 'FontSize', 14)
    end
end

%% Mean plots
if plotmeans
    % Control task
    figure('WindowStyle', 'docked')
    ylims = [300 1000];
    elines = {' sk', ' dk', ' ok'};
    lines = {':sk', '--dk', '-ok'};
    plines = {' sr', ' dr', ' or'};
    subplot(2,3,1)
    items = [1 3 5];
    for i = 1:numel(items)
        iset = mm{1}(mm{1}(:,2) == items(i), :);
        sset = sm{1}(mm{1}(:,2) == items(i), :);
%         pset = pm{1}(pm{1}(:,2) == items(i), :);
        errorbar(1:6, iset(:,3), sset(:,3)./sqrt(8), elines{i}, 'LineWidth', 2, 'MarkerSize', 10);
        hold on
        h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
%         plot(1:6, pset(:,3), plines{i}, 'LineWidth', 2)
        
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('Mean RT', 'FontSize', 14)
    title('Control')
    legend(h, 'Far', 'Mid', 'Near')
    
    % Correlated task
    lines = {':sk', '--dk', '-ok'};
    subplot(2,3,2)
    items = [1 3 5];
    for i = 1:numel(items)
        iset = mm{2}(mm{2}(:,2) == items(i), :);
        sset = sm{2}(mm{2}(:,2) == items(i), :);
%         pset = pm{2}(pm{2}(:,2) == items(i), :);
        errorbar(1:6, iset(:,3), sset(:,3)./sqrt(8), elines{i}, 'LineWidth', 2, 'MarkerSize', 10);
        hold on
        h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
%         plot(1:6, pset(:,3), plines{i}, 'LineWidth', 2)
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('Mean RT', 'FontSize', 14)
    title('Correlated')
    legend(h, 'Far', 'Mid', 'Near')
    
    % Filtering task
    elines{1} = {' sk', ' dk', ' ok'};
    elines{2} = {' sk', ' dk', ' ok'};
    lines{1} = {':sk', '--dk', '-ok'};
    lines{2} = {':sk', '--dk', '-ok'};
    plines{1} = {' sr', ' dr', ' or'};
    plines{2} = {' sr', ' dr', ' or'};
    
    offsets = [-.15, .15];
    subplot(2,3,3)
    items = [1 3 5];
    for i = 1:numel(items)
        if i == 3
            item = items(i);
            sets = [1:2:12; 2:2:12];
            for j = 1:2
                iset = mm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
                sset = sm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
                %             pset = pm{3}(pm{3}(:,2) == items(i) & ismember(pq{3}(:,1), sets(j,:)), :);
                
                errorbar((1:6)+offsets(j), iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                hold on
                if j == 1
                    h(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                else
                    h2(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                end
            end
        else
            item = items(i);
            sets = [1:2:12];
            for j = 1
                iset = mm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
                sset = sm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
                %             pset = pm{3}(pm{3}(:,2) == items(i) & ismember(pq{3}(:,1), sets(j,:)), :);
                
                errorbar((1:6)+offsets(j), iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                hold on
                if j == 1
                    h(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                else
                    plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                end
            end
        end
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('Mean RT', 'FontSize', 14)
    title('Filtering')
    legend([h, h2(3)], 'Far', 'Mid', 'Near - Same IDV', 'Near - Different IDV')
end

%% Plot accuracy
if plotaccuracy
    % Control task
%     figure('WindowStyle', 'docked')
    ylims = [0 .5];
    elines = {' sk', ' dk', ' ok'};
    lines = {':sk', '--dk', '-ok'};
%     plines = {' sr', ' xr', ' or'};
    
    subplot(2,3,3+1)
    items = [1 3 5];
    for i = 1:numel(items)
        iset = ma{1}(ma{1}(:,2) == items(i), :);
        sset = sa{1}(ma{1}(:,2) == items(i), :);
%         pset = pm{1}(pm{1}(:,2) == items(i), :);
        errorbar(1:6, iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{i}, 'LineWidth', 2, 'MarkerSize', 10);
        hold on
        h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
%         plot(1:6, pset(:,3), plines{i}, 'LineWidth', 2)
        
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('p(Error)', 'FontSize', 14)
    title('Control')
    legend(h, 'Far', 'Mid', 'Near')
    
    % Correlated task
    lines = {':sk', '--dk', '-ok'};
    subplot(2,3,3+2)
    items = [1 3 5];
    for i = 1:numel(items)
        iset = ma{2}(ma{2}(:,2) == items(i), :);
        sset = sa{2}(ma{2}(:,2) == items(i), :);
%         pset = pm{2}(pm{2}(:,2) == items(i), :);
        errorbar(1:6, iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{i}, 'LineWidth', 2, 'MarkerSize', 10);
        hold on
        h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
%         plot(1:6, pset(:,3), plines{i}, 'LineWidth', 2)
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('p(Error)', 'FontSize', 14)
    title('Correlated')
    legend(h, 'Far', 'Mid', 'Near')
    
 % Filtering task
    elines{1} = {' sk', ' dk', ' ok'};
    elines{2} = {' sk', ' dk', ' ok'};
    lines{1} = {':sk', '--dk', '-ok'};
    lines{2} = {':sk', '--dk', '-ok'};
    plines{1} = {' sr', ' dr', ' or'};
    plines{2} = {' sr', ' dr', ' or'};
    
    offsets = [-.15, .15];
    subplot(2,3,3+3)
    items = [1 3 5];
    for i = 1:numel(items)
        if i == 3
            item = items(i);
            sets = [1:2:12; 2:2:12];
            for j = 1:2
                iset = ma{3}(ma{3}(:,2) == items(i) & ismember(ma{3}(:,1), sets(j,:)), :);
                sset = sa{3}(ma{3}(:,2) == items(i) & ismember(ma{3}(:,1), sets(j,:)), :);
                %             pset = pm{3}(pm{3}(:,2) == items(i) & ismember(pq{3}(:,1), sets(j,:)), :);
                
                errorbar((1:6)+offsets(j), iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                hold on
                if j == 1
                    h(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                else
                    h2(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                end
            end
        else
            item = items(i);
            sets = [1:2:12];
            for j = 1
                iset = ma{3}(ma{3}(:,2) == items(i) & ismember(ma{3}(:,1), sets(j,:)), :);
                sset = sa{3}(ma{3}(:,2) == items(i) & ismember(ma{3}(:,1), sets(j,:)), :);
                %             pset = pm{3}(pm{3}(:,2) == items(i) & ismember(pq{3}(:,1), sets(j,:)), :);
                
                errorbar((1:6)+offsets(j), iset(:,3), sset(:,3)./sqrt(numel(subIdxs)), elines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                hold on
                if j == 1
                    h(i) = plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                else
                    plot((1:6)+offsets(j), iset(:,3), lines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
                    %                 plot((1:6)+offsets(j), pset(:,3), plines{j}{i}, 'LineWidth', 2)
                end
            end
        end
    end
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('p(Error)', 'FontSize', 14)
    title('Filtering')
    legend([h, h2(3)], 'Far', 'Mid', 'Near - Same IDV', 'Near - Different IDV')
end
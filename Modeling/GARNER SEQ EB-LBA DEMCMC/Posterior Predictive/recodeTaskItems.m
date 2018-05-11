% Recode control2, correlated2, and the even numbered filtration items 
function taskDataMatrix = recodeTaskItems(taskDataMatrix, cols)

itemSets = {1:2:12; 2:2:12; [1 3 5 8 10 12]; [2 4 6 7 9 11]; 1:12}; % Sets of items for each condition
task = unique(taskDataMatrix(:,strcmp(cols, 'task')));

switch task
    case {1, 3}
        % Do nothing
    case {2, 4}
        % Recode 
        items = itemSets{task};
        recode = itemSets{task - 1};
        for i = 1:numel(items)
            taskDataMatrix(taskDataMatrix(:,strcmp(cols, 'item')) == items(i), strcmp(cols,'item')) = recode(i);
            taskDataMatrix(taskDataMatrix(:,strcmp(cols, 'prev')) == items(i), strcmp(cols,'prev')) = recode(i);
        end
    case 5
        % Recode irrelevant dimension change
        target = [2, 4, 6];
        oldmap = [2 4 6 8 10 12 1 3 5 7 9 11];
        newmap = circshift(oldmap, [0, 6]); 
        f = @(x)(find(x == oldmap));
        
        % Recode previous item if current item is even
        taskDataMatrix(ismember(taskDataMatrix(:,strcmp(cols, 'item')), target), strcmp(cols, 'prev')) =...
            newmap(arrayfun(f, taskDataMatrix(ismember(taskDataMatrix(:,strcmp(cols, 'item')), target), strcmp(cols, 'prev'))))';
        
        % Recode current item if current item is even
        taskDataMatrix(ismember(taskDataMatrix(:,strcmp(cols, 'item')), target), strcmp(cols, 'item')) =...
            taskDataMatrix(ismember(taskDataMatrix(:,strcmp(cols, 'item')), target), strcmp(cols, 'item')) - 1;
end    
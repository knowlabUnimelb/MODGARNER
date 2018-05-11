% This function will set up the previous trial item code
% Also, for plotting, this will flip the indexes so they refer only to
% near, middle, and far rather than A and B
function [taskDataMatrix, cols] = codePreviousTrial(taskDataMatrix, cols)

itemSets = {1:2:12; 2:2:12; [1 3 5 8 10 12]; [2 4 6 7 9 11]; 1:12}; % Sets of items for each condition

% prev = taskDataMatrix(taskDataMatrix(:,strcmp(cols, 'trl')) ~= 120, strcmp(cols, 'item'));
% curr = taskDataMatrix(taskDataMatrix(:,strcmp(cols, 'trl')) ~= 1, :);
prev = taskDataMatrix(1:end-1, strcmp(cols, 'item'));
curr = taskDataMatrix(2:end, :);
task = unique(taskDataMatrix(:,strcmp(cols, 'task')));

%% 
switch task
    case {1, 2, 3, 4} % Control or correlated conditions
        rev = flip(itemSets{task});        
    case 5 % Filtering conditions
        space = reshape(itemSets{task}, 2, 6);
        temp = fliplr(space);
        rev = temp(:);
end

% Need to flip B items so they are coded as A items
% Prev needs to be coded as far_same, mid_same, near_same, near_opp, mid_opp, far_opp
% Curr needs to be codes as far, mid, same
rev_items = itemSets{task}(itemSets{task} >= 7);
findfun = @(x)find(x==itemSets{task}); % Function to get the indices in itemSets{task}
findrev = @(x)find(x==rev); % Function to get the indices in rev
for i = 1:numel(rev_items) % Cycle through category B items
    % Reverse previous item indexes for category B items
    prev(curr(:,strcmp(cols, 'item')) == rev_items(i)) =... % For any of the previous items when the current item is from category B
        rev(arrayfun(findfun, prev(curr(:,strcmp(cols, 'item')) == rev_items(i)))); % Recode the previous item
    
    % Reverse category B items
    curr(curr(:,strcmp(cols, 'item')) == rev_items(i), strcmp(cols, 'item')) =...
        itemSets{task}(arrayfun(findrev, curr(curr(:,strcmp(cols, 'item')) == rev_items(i), strcmp(cols, 'item'))));
end    

% Add columns back into taskDataMatrix
taskDataMatrix = [curr, prev];
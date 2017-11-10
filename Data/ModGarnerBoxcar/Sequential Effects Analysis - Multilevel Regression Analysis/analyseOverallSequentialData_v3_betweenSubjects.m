% v2 - test against mean of other items 
% v3 - test against mean of other items excluding repetitions

clear
clc
close all

whichCondition = 'Sat'; % {'Lin', 'Sat'}
whichtask = 1:4; % 1-4

%% Test data
cols = {'sub', 'task', 'prev', 'item', 'rt'};
mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\MODIFIED GARNER BOXCAR\';
dataloc = fullfile(mainfolder, 'Sequential Effects Analysis - Organize Data');
data =  dlmread(fullfile(dataloc, sprintf('test_boxcar_seq%sCorrectRT.dat', whichCondition)));
taskNames = {'Control', 'Correlated', 'Filtration', 'FiltrationIDV'};

savename = sprintf('Garner%s%s', whichCondition, taskNames{whichtask});

sub = data(:,strcmp(cols, 'sub'));
ns = numel(unique(sub));

% Recode items so that the filtration other irrelevant dimension equals
% task 4 and that items are numbered 1-3 and previous items are numbered
% 1-6
temp = data(:,strcmp(cols, 'task'));
temp(mod(data(:,strcmp(cols, 'prev')), 2) == 0) = 4;
data(:,strcmp(cols, 'task')) = temp;
data(data(:,strcmp(cols, 'task')) == 4, strcmp(cols, 'prev')) =...
    data(data(:,strcmp(cols, 'task')) == 4, strcmp(cols, 'prev')) - 1;

% data(:,strcmp(cols, 'prev')) = floor(data(:,strcmp(cols, 'prev'))./2) + 1;
% data(:,strcmp(cols, 'item')) = floor(data(:,strcmp(cols, 'item'))./2) + 1;

sub = [1 2 3 4];

data = data(ismember(data(:,strcmp(cols, 'sub')), sub) & ismember(data(:,strcmp(cols, 'task')), whichtask), :);

pi = unique(data(:,2:4), 'rows');
ui = unique(data(:,strcmp(cols, 'item')));
up = unique(data(:,strcmp(cols, 'prev')));
ut = unique(data(:,strcmp(cols, 'task')));

temp = data;
tcnt = 1;
for tidx = 1:numel(ut)
    pcnt = 1;
    for pidx = 1:numel(up)
        icnt = 1;
        for iidx = 1:numel(ui);
            temp(data(:,strcmp(cols, 'task')) == ut(tidx) &...
                data(:,strcmp(cols, 'prev')) == up(pidx) &...
                data(:,strcmp(cols, 'item')) == ui(iidx), strcmp(cols, 'task')) = tcnt;
            temp(data(:,strcmp(cols, 'task')) == ut(tidx) &...
                data(:,strcmp(cols, 'prev')) == up(pidx) &...
                data(:,strcmp(cols, 'item')) == ui(iidx), strcmp(cols, 'prev')) = pcnt;
            temp(data(:,strcmp(cols, 'task')) == ut(tidx) &...
                data(:,strcmp(cols, 'prev')) == up(pidx) &...
                data(:,strcmp(cols, 'item')) == ui(iidx), strcmp(cols, 'item')) = icnt;
            icnt = icnt + 1;
        end
        pcnt = pcnt + 1;
    end
    tcnt = tcnt + 1;
end
data = temp;

%%
subs = data(:,strcmp(cols, 'sub'));
ns = numel(unique(subs));

task = data(:, strcmp(cols, 'task'));
nt = numel(unique(task));

prev = data(:,strcmp(cols, 'prev'));
np = numel(unique(prev));

item = data(:,strcmp(cols, 'item'));
ni = numel(unique(item));

rt = data(:,strcmp(cols, 'rt'));
nd = numel(rt);

d = [subs task item prev rt];

%% Stats

testItems = [3; 1; 4];
conditions = {'Control', 'Correlated', 'Filtration'};
testNames = {'Boundary vs Mean RT', 'Far vs Mean RT', 'Adjacent vs Mean RT'};
mldata = [];
ns = [];
for i = 1:3 % Task
    for j = 1:3 % Item
        d1 = d( d(:,2) == i & d(:,3) == 3 & d(:,4) == testItems(j,1), [1 5]);
        if j == 1
            d2 = d( d(:,2) == i & d(:,3) == 3 & d(:,4) ~= testItems(j,1), [1 5]);
        else
            d2 = d( d(:,2) == i & d(:,3) == 3 & d(:,4) ~= testItems(j,1) & d(:,4) ~= 3, [1 5]);
        end
        
       
        mldata = [mldata; 
         d1(:,1), i * ones(size(d1,1),1), j * ones(size(d1, 1), 1), 1 * ones(size(d1, 1), 1), d1(:,2);
         d2(:,1), i * ones(size(d2,1),1), j * ones(size(d2, 1), 1), 0 * ones(size(d2, 1), 1), d2(:,2)];
             ns = [ns; size(d1,1) size(d2,1)];
    end
end

% Test irrelevant dimension
d1 = d( d(:,2) == 3 & d(:,3) == 3 & d(:,4) == 3, [1 5]); % Filtration, near current, near prev
d2 = d( d(:,2) == 4 & d(:,3) == 3 & d(:,4) == 3, [1 5]); % Filtration IDV, near current, near prev
d3 = d( d(:,2) == 3 & d(:,3) == 3 & d(:,4) == 1, [1 5]); % Filtration, near current, far prev
d4 = d( d(:,2) == 4 & d(:,3) == 3 & d(:,4) == 1, [1 5]); % FiltrationIDV, near current, far prev
d5 = d( d(:,2) == 3 & d(:,3) == 3 & d(:,4) == 4, [1 5]); % Filtration, near current, adj prev
d6 = d( d(:,2) == 4 & d(:,3) == 3 & d(:,4) == 4, [1 5]); % FiltrationIDV, near current, far prev

mldata2 = [d1(:,1), 3 * ones(size(d1,1),1), 1 * ones(size(d1, 1), 1), d1(:,2); % filt, near prev
           d2(:,1), 4 * ones(size(d2,1),1), 1 * ones(size(d2, 1), 1), d2(:,2); % filtIDV, near prev
           d3(:,1), 3 * ones(size(d3,1),1), 2 * ones(size(d3, 1), 1), d3(:,2); % filt, far prev
           d4(:,1), 4 * ones(size(d4,1),1), 2 * ones(size(d4, 1), 1), d4(:,2); % filtIDV, far prev
           d5(:,1), 3 * ones(size(d5,1),1), 3 * ones(size(d5, 1), 1), d5(:,2); % filt, adj prev
           d6(:,1), 4 * ones(size(d6,1),1), 3 * ones(size(d6, 1), 1), d6(:,2)]; % filtIDV, adj prev;


if strcmp(whichCondition, 'Lin')
    dlmwrite('mldata_boxcar_rt_lin.dat', [1*ones(size(mldata, 1), 1), mldata], '\t')
    dlmwrite('mldata_boxcar_rt_idv_lin.dat', [1*ones(size(mldata2, 1), 1), mldata2], '\t')
else
    dlmwrite('mldata_boxcar_rt_sat.dat', [2*ones(size(mldata, 1), 1), mldata], '\t')
    dlmwrite('mldata_boxcar_rt_idv_sat.dat', [2*ones(size(mldata2, 1), 1), mldata2], '\t')
end
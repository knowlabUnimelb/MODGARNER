clear all
clc

subject = 204; % Odd is brightness/even is saturation 101:104, 201:204
condition = 2;
sessions = 2:6; 
lastsession = 6;

%SEPARABLE
cols = {'sub', 'con', 'side', 'sess', 'task', 'prac', 'trl', 'itm', 'leftBri', 'leftSat', 'rightBri', 'rightSat', 'corrResp', 'resp', 'corr', 'rt', 'exc'};

%INTEGRAL
%cols = {'sub', 'con', 'sess', 'task', 'prac', 'trl', 'itm', 'sat', 'bri', 'corrResp', 'resp', 'corr', 'rt', 'exc'};

mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\MODIFIED GARNER SEPARABLE';
dataloc = fullfile(mainfolder, 'rawdata');

data = [];
for i = 1:numel(sessions)
    %SEPARABLE
    data = [data; dlmread(fullfile(dataloc, sprintf('GARNER_SEPARABLE_%03d_%02d_s%02d.dat', subject, condition, sessions(i))))];
end
data(data(:,strcmp(cols, 'task')) == 6, strcmp(cols, 'task')) = 5; 

%% Setup output label
if ismember(subject, [101, 102, 103, 104])
    outfileidx = sprintf('B%d', find(subject == [101, 102, 103, 104]));
elseif ismember(subject, [201, 202, 203, 204])
    outfileidx = sprintf('S%d', find(subject == [201, 202, 203, 204]));
else 
    error('Bad subject number')    
end

%% Remove outliers
% Remove rts < 200 & > 3000
% data(data(:,strcmp(cols, 'corr')) == 9, :) = [];
data(data(:, strcmp(cols, 'prac'))==1 & data(:,strcmp(cols, 'trl')) <= 24, :) = []; % Remove practice
alldata = data;

boundaryitems = [5 6 7 8];
badrts = data(data(:,strcmp(cols, 'rt')) < .2 | data(:,strcmp(cols, 'rt')) > 3, :);
nbadrts = size(badrts, 1);

%150626 Don't remove outliers but create a column to tag them
%data(data(:,strcmp(cols, 'rt')) < 200 | data(:,strcmp(cols, 'rt')) > 3000, :) = [];
rtOutliers = find(data(:,strcmp(cols, 'rt')) < .2 | data(:,strcmp(cols,'rt')) > 3);

timeouts = data(data(:,strcmp(cols, 'corr')) == 9, :);
timeOutliers = find(data(:,strcmp(cols, 'corr')) == 9);

exclude = zeros(size(data, 1),1);
data = [data, exclude];
data(rtOutliers, strcmp(cols, 'exc')) = 1;
data(timeOutliers, strcmp(cols, 'exc')) = 1;

cdata = []; acc = nan(12, 5); mrt = nan(12,5); emrt = nan(12, 5); smrt = nan(12, 5);
for i = 1:numel(unique(data(:, strcmp(cols, 'task'))))
   tdata = data(data(:,strcmp(cols, 'task')) == i, :);
   items{i} = unique(tdata(:,strcmp(cols, 'itm')));
   for j = 1:numel(items{i})
       acc(items{i}(j),i) = nanmean(tdata(tdata(:, strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'exc')) == 0, strcmp(cols, 'corr')));
       mrt(items{i}(j),i) = nanmean(tdata(tdata(:, strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 1 & tdata(:,strcmp(cols, 'exc')) == 0, strcmp(cols, 'rt')));
       smrt(items{i}(j),i) = nanstd(tdata(tdata(:, strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 1 & tdata(:,strcmp(cols, 'exc')) == 0, strcmp(cols, 'rt')));
       
       outliers{items{i}(j),i} = tdata(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 1 & tdata(:,strcmp(cols, 'exc')) == 0 & tdata(:,strcmp(cols, 'rt')) > 3 * smrt(items{i}(j),i) + mrt(items{i}(j),i), :);
       outlierIdx{items{i}(j),i} = find(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 1 & tdata(:,strcmp(cols, 'exc')) == 0 & tdata(:,strcmp(cols, 'rt')) > 3 * smrt(items{i}(j),i) + mrt(items{i}(j),i));
       noutliers(items{i}(j),i) = size(outliers{items{i}(j),i}, 1);
       
       %150626 Don't remove outliers but create a column to tag them
       %tdata(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 1 & tdata(:,strcmp(cols, 'rt')) > (3 * smrt(items{i}(j),i) + mrt(items{i}(j),i)), :) = [];
       tdata(outlierIdx{items{i}(j),i}, strcmp(cols, 'exc')) = 1;
       
       boundaryMRT{i} = nanmean(mrt(items{i}(ismember(items{i}, boundaryitems)), i));
%        boundaryERR{i} = 1 - mean(acc(items{i}(ismember(items{i}, boundaryitems)), i));
       boundaryERRa{i} = 1 - nanmean(acc(items{i}(ismember(items{i}, boundaryitems(1:2))), i));
       boundaryERRb{i} = 1 - nanmean(acc(items{i}(ismember(items{i}, boundaryitems(3:4))), i));
       
       emrt(items{i}(j),i) = nanmean(tdata(tdata(:, strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 0 & tdata(:,strcmp(cols, 'exc')) == 0, strcmp(cols, 'rt')));
       esmrt(items{i}(j),i) = nanstd(tdata(tdata(:, strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 0 & tdata(:,strcmp(cols, 'exc')) == 0, strcmp(cols, 'rt')));
       eoutliers{items{i}(j),i} = tdata(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 0 & tdata(:,strcmp(cols, 'exc')) == 0 & tdata(:,strcmp(cols, 'rt')) > 3 * esmrt(items{i}(j),i) + emrt(items{i}(j),i), :);
       eoutlierIdx{items{i}(j),i} = find(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 0 & tdata(:,strcmp(cols, 'exc')) == 0 & tdata(:,strcmp(cols, 'rt')) > 3 * esmrt(items{i}(j),i) + emrt(items{i}(j),i));
       enoutliers(items{i}(j),i) = size(eoutliers{items{i}(j),i}, 1);
       
       %150626 Don't remove outliers but create a column to tag them
       %tdata(tdata(:,strcmp(cols, 'itm')) == items{i}(j) & tdata(:,strcmp(cols, 'corr')) == 0 & tdata(:,strcmp(cols, 'rt')) > (3 * esmrt(items{i}(j),i) + emrt(items{i}(j),i)), :) = [];
       tdata(eoutlierIdx{items{i}(j),i}, strcmp(cols, 'exc')) = 1;
   end
   cdata = [cdata; tdata];
end
percentRemoved = (size(data, 1) - size(cdata, 1))./size(data, 1) * 100; 

%% Convert to RT
cdata(:,strcmp(cols, 'rt')) = cdata(:,strcmp(cols, 'rt')) * 1000; 

%%
% fdata = sortrows(cdata(:, mstrfind(cols, {'sub', 'task',  'itm', 'corr', 'rt'})), [1 2 3 4]); 
fdata0 = cdata(:, mstrfind(cols, {'sub', 'task', 'sess', 'trl', 'itm', 'corr', 'rt', 'exc'})); 
fdata = [];
sesstaskorder = data(data(:,strcmp(cols, 'trl')) == 1,4:5); %for integral, change to 3:4
for i = 2:lastsession
    idx(i-1) = find(sesstaskorder(:,1) == i & sesstaskorder(:,2) == 5, 1, 'last');
end
sesstaskorder(idx,:) = [];

for i = 1:size(sesstaskorder,1)
        fdata = [fdata; fdata0(fdata0(:,2) == sesstaskorder(i,2) & fdata0(:,3) == sesstaskorder(i,1),:)];
end

%% Write data to new file
fid = fopen(sprintf('Garner_SEP_sub%s_ordered_tagged_msec.dat', outfileidx), 'wt');
% fprintf(fid, '%3s\t %4s\t %4s\t %3s\t %4s\t %3s\t %7s\t %3s\n', 'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc');
fprintf(fid, '%4d\n', size(fdata,1));
for i = 1:size(fdata, 1)
    fprintf(fid, '%3d\t %4d\t %4d\t %3d\t %4d\t %3d\t %7.2f\t %3d\n', fdata(i,1), fdata(i,2), fdata(i,3), fdata(i,4), fdata(i,5), fdata(i,6), fdata(i,7), fdata(i,8));
end
fclose(fid);

disp('done')
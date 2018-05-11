clear all
clc

% OLD
% load testdata
% 
% n = 5; % size(simrt, 2);
dataprefix = 'Garner_SIM2_subS%02d_ordered_tagged.dat';
% 
% for i = 1:n
%    simdata = data;
%    simdata(:,strcmp(cols, 'rts')) = simrt(:,i) * 1000;
%    simdata(:,strcmp(cols, 'acc')) = simacc(:,i);
%    simdata(:,strcmp(cols, 'sub')) = 100+i;
%    dlmwrite(sprintf(dataprefix, i), 3600)
%    dlmwrite(sprintf(dataprefix, i), simdata(:,1:8), '-append', 'delimiter', '\t')
% end


% NEW    
mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\GarnerDEMCMC_recode';
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc', 'igp', 'cat', 'opp'}; % Data columns

%            w   wc     c    cf    b1    b2    b3    b4    b5    wp    bh    sp     A     B    vs    t0
parms = [0.35, 0.50, 0.15, 0.10, 0.50, 0.50, 0.50, 0.50, 0.50, 0.80, 1.50, 0.14, 0.27, 0.27, 0.20, 0.10];
stimloc = load(fullfile(mainfolder, 'MDS', sprintf('mds_Separable.dat')));
data = dlmread(fullfile(mainfolder, 'Data', 'Garner_BOXCAR_subL1_ordered_tagged.dat'));
data(1,:) = [];
data(:,strcmp(cols, 'igp'))  = data(:, strcmp(cols, 'trl')) == 1; % Equals 1 if it's the start of a new block
data(:,strcmp(cols, 'cat')) = ones(size(data, 1), 1);  % Category for each item
data(data(:,strcmp(cols, 'item')) >= 7, strcmp(cols, 'cat')) = 2;
data(:,strcmp(cols, 'opp')) = 2 - (data(:,strcmp(cols, 'cat')) - 1);             % Opposite category for each item
  
for i = 1:5
    sdata = data; 
    [rt, acc] = seqEBLBA(parms, stimloc, sdata, cols, true);
    sdata(:,strcmp(cols, 'rts')) = rt * 1000;
    sdata(:,strcmp(cols, 'acc')) = acc;
    sdata(:, strcmp(cols, 'exc')) = zeros(numel(rt, 1), 1);
    sdata(rt < .2 | rt > 3, strcmp(cols, 'exc')) = 1;
    sdata(:,strcmp(cols, 'sub')) = 200 + i;
    dlmwrite(sprintf(dataprefix, i), 3600)
    dlmwrite(sprintf(dataprefix, i), sdata(:,1:8), '-append', 'delimiter', '\t')
end
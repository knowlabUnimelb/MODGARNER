clear all
clc
close all

mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\MODIFIED GARNER SEPABLE';
datafolder = fullfile(mainfolder, 'data');
filenameStr = 'Garner_SEP_sub%s_ordered_tagged_msec.dat';
subjects = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'};
% filenameStr = 'Garner_BOXCAR_sub%s_ordered_tagged_msec.dat';
% subjects = {'L1', 'L2', 'L3', 'L4', 'S1', 'S2', 'S3', 'S4'};
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rt', 'exc'};

cnt = 1;
for i = 1:numel(subjects)
   % Load subject data
   data{cnt} = dlmread(fullfile(datafolder, sprintf(filenameStr, subjects{i})));
   data{cnt}(1,:) = []; % Remove first line which indicates the number of rows in the data file
   
   % Recode tasks
   tasks = data{cnt}(:,strcmp(cols, 'task'));
   data{cnt}(ismember(tasks, 1:2), strcmp(cols, 'task')) = 1;
   data{cnt}(ismember(tasks, 3:4), strcmp(cols, 'task')) = 2;
   data{cnt}(ismember(tasks, 5), strcmp(cols, 'task'))   = 3;
   
   % Recode items
   items = data{cnt}(:,strcmp(cols, 'item'));
   data{cnt}(ismember(items, [5 6 7 8]), strcmp(cols, 'item')) = 1;
   data{cnt}(ismember(items, [3 4 9 10]), strcmp(cols, 'item')) = 2;
   data{cnt}(ismember(items, [1 2 11 12]), strcmp(cols, 'item'))= 3;
   cnt = cnt + 1;
end

% Remove outliers
data = cellfun(@(x)(x(~x(:,strcmp(cols, 'exc')), :)), data, 'UniformOutput', false);

% Overall Averages
meanOverallAccuracy  = cellfun(@(x)(aggregate(x                                  ,...
    mstrfind(cols, {'task'}), mstrfind(cols, {'acc'}), @mean, 1)), data, 'UniformOutput', false);
meanOverallCorrectRT = cellfun(@(x)(aggregate(x(x(:,strcmp(cols, 'acc')) == 1, :),...
    mstrfind(cols, {'task'}), mstrfind(cols, {'rt'}), @mean, 1)), data, 'UniformOutput', false);

stdOverallAccuracy  = cellfun(@(x)(aggregate(x                                  ,...
    mstrfind(cols, {'task'}), mstrfind(cols, {'acc'}), @std, 1)), data, 'UniformOutput', false);
stdOverallCorrectRT = cellfun(@(x)(aggregate(x(x(:,strcmp(cols, 'acc')) == 1, :),...
    mstrfind(cols, {'task'}), mstrfind(cols, {'rt'}), @std, 1)), data, 'UniformOutput', false);

nOverallAccuracy  = cellfun(@(x)(aggregate(x                                  ,...
    mstrfind(cols, {'task'}), mstrfind(cols, {'acc'}), @count, 1)), data, 'UniformOutput', false);
nOverallCorrectRT = cellfun(@(x)(aggregate(x(x(:,strcmp(cols, 'acc')) == 1, :),...
    mstrfind(cols, {'task'}), mstrfind(cols, {'rt'}), @count, 1)), data, 'UniformOutput', false);

%% JAGS analysis

%% MCMC Settings for JAGS
nchains  = 2;     % How Many Chains?
nburnin  = 500;   % How Many Burn-in Samples?
nsamples = 5000;  % How Many Recorded Samples?

%% Assign Matlab Variables to the Observed JAGS Nodes
nConditions = size(cell2mat(meanOverallCorrectRT),1);
nSubjects = size(cell2mat(meanOverallCorrectRT),2);

deviation = cell2mat(meanOverallCorrectRT) - ([1 1 1]' * mean(cell2mat(meanOverallCorrectRT)));

datastruct = struct('subD', deviation', 'c', nConditions, 'n', nSubjects);

%% Initialize Unobserved Variables
% Starting points for the unobserved variables in JAGS
for i=1:nchains
    S.mu = zeros(nConditions, 1);             % An intial Value for the Gaussian likelihood mean
    S.sigma = ones(nConditions,1); % Intial values For All The precisions
    M = 0;
    L = .001;
    init0(i) = S;
end

%% Run JAGS
jagsModelFileName = 'overallRTmodel.txt';
whichParmsToMonitor =  {'mu','sigma', 'lambda', 'predD'};

% Pass all information to JAGS
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, jagsModelFileName), ...
    init0, ...
    'doparallel' , false, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', 10, ...
    'monitorparams', whichParmsToMonitor, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0 , ...
    'workingdir' , 'tmpjags' );

%% Plot subject samples
for idx = 1:nSubjects; 
    d = deviation(:,idx);
    plotChains(squeeze(samples.predD(:,:,idx,:)),...
        {sprintf('con = %3.2f', d(1)),...
        sprintf('corr = %3.2f', d(2)),...
        sprintf('filt = %3.2f', d(3))});
end

%% Plot condition means
plotChains(samples.mu, {'con', 'corr', 'filt'})
plotChains(samples.sigma, {'con', 'corr', 'filt'})

%%
figure
muSamp = [squeeze(samples.mu(1,:,:)); squeeze(samples.mu(2,:,:))];
bandwidth = max(getbandwidth(muSamp(:)), .01);
[h, L, MX, MED, bw] = violin(muSamp,...
     'facecolor', [1 1 1], 'edgecolor', 'k',...
     'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
 hold on
 line([.5 3.5], [0 0], 'LineStyle', '--', 'Color', 'k')
 set(gca,'XTick', 1:3, 'XTickLabel', {'Control', 'Correlated', 'Filtering'})
 xlabel('Condition')
 ylabel('Deviation from Overall Mean')
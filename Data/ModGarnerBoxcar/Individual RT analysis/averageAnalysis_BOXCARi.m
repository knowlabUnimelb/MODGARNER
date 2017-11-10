clear all
clc
close all

mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\MODIFIED GARNER BOXCAR';
datafolder = fullfile(mainfolder, 'data');
filenameStr = 'Garner_BOXCAR_sub%s_ordered_tagged_msec.dat';
subjects = {'L1', 'L2', 'L3', 'L4', 'S1', 'S2', 'S3', 'S4'};
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
jagsdata = cell(1,8);
for subIdx = 1:8
    jagsdata{subIdx} = [];
    for taskIdx = 1:3
        for itemIdx = 1:3
            jagsdata{subIdx} = [jagsdata{subIdx};
                data{subIdx}(data{subIdx}(:,strcmp(cols, 'task')) == taskIdx &...
                   data{subIdx}(:,strcmp(cols, 'item')) == itemIdx &...
                    data{subIdx}(:,strcmp(cols, 'acc')) == 1, strcmp(cols, 'task')),...
                data{subIdx}(data{subIdx}(:,strcmp(cols, 'task')) == taskIdx &...
                    data{subIdx}(:,strcmp(cols, 'item')) == itemIdx &...
                    data{subIdx}(:,strcmp(cols, 'acc')) == 1, strcmp(cols, 'item')),...
                log(data{subIdx}(data{subIdx}(:,strcmp(cols, 'task')) == taskIdx &...
                    data{subIdx}(:,strcmp(cols, 'item')) == itemIdx &...
                    data{subIdx}(:,strcmp(cols, 'acc')) == 1, strcmp(cols, 'rt')))];
        end
    end
end

meanOverallLogRT = cell2mat(cellfun(@(x)(aggregate(x, 1:2, 3, @mean, 1)), jagsdata, 'UniformOutput', false));

%% MCMC Settings for JAGS
nchains  = 2;     % How Many Chains?
nburnin  = 500;   % How Many Burn-in Samples?
nsamples = 5000;  % How Many Recorded Samples?

for fitsub = 1:8;

task = jagsdata{fitsub}(:,1);
item = jagsdata{fitsub}(:,2);
logrt = jagsdata{fitsub}(:,3);

%% Assign Matlab Variables to the Observed JAGS Nodes
datastruct = struct('logrt', logrt,...
                    'task', task,...
                    'item', item,...
                    'nRT', size(jagsdata{fitsub}, 1),...
                    'nTask', 3,...
                    'nItem', 3);

%% Initialize Unobserved Variables
% Starting points for the unobserved variables in JAGS
for i=1:nchains
    S.mu = 6 * ones(3, 3);             % An intial Value for the Gaussian likelihood mean
    S.sigma = ones(3, 3); % Intial values For All The precisions
    M = 6;
    L = .001;
    init0(i) = S;
end

%% Run JAGS
jagsModelFileName = 'individualLogRTmodel.txt';
whichParmsToMonitor =  {'mu','sigma', 'lambda', 'predrt'};

% Pass all information to JAGS
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, jagsModelFileName), ...
    init0, ...
    'doparallel' , true, ...
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
% for idx = 1:nSubjects; 
%     d = deviation(:,idx);
%     plotChains(squeeze(samples.predD(:,:,idx,:)),...
%         {sprintf('con = %3.2f', d(1)),...
%         sprintf('corr = %3.2f', d(2)),...
%         sprintf('filt = %3.2f', d(3))});
% end

% %% Plot condition means
plotChains(samples.mu, {'con, near', 'con, mid', 'con, far', 'corr, near', 'corr, mid', 'corr, far', 'filt, near', 'filt, mid', 'filt, far'})
plotChains(samples.sigma, {'con, near', 'con, mid', 'con, far', 'corr, near', 'corr, mid', 'corr, far', 'filt, near', 'filt, mid', 'filt, far'})
% 
%%
% predrt = [squeeze(samples.predrt(1,:,:)); squeeze(samples.predrt(2,:,:))]';
% 

mfun    = @(mu, sigma)(exp(sigma.^2./2 + mu));
figure
adj = [-.5 0 .5];
colors = [1 0 0; 0 1 0; 0 0 1];
for i = 1:3
    pRT = nan(1000, 3);
    for j = 1:3
% %         prt = predrt(task == i & item == j, :);
%         prt = exp(logrt(task == i & item == j));
        postmu = squeeze(samples.mu(:,:,i,j));
        postsigma = squeeze(samples.sigma(:,:,i,j));
        tempmu = postmu(:);
        tempsig = postsigma(:);
        postm = mfun(tempmu, tempsig);
        pRT(:,j) = randsample(postm, 1000, false); % sample without replacement

        bandwidth = max(getbandwidth(pRT(:)), .01);
    end
%         
    [h, L, MX, MED, bw] = violin(pRT,...
        'x', (2:2:6) + adj(i),...
        'facecolor', colors(i,:), 'edgecolor', 'k',...
        'facealpha', .5,...
        'mc', '', 'medc', '', 'bw', bandwidth);
    hold on
    plot((2:2:6) + adj(i), mean(pRT), '-k', 'LineWidth', 1)
    plot((2:2:6) + adj(i), mean(pRT), ' ok')
    vh(i) = h(1);
    
    means(i,:) = mean(pRT);
    cis(:,:,i) = prctile(pRT, [.5 99.5]);
end
set(gca,'XTick', 2:2:6, 'XTickLabel', {'Near', 'Mid', 'Far'}, 'XLim', [0 8], 'YLim', [300 1000])
xlabel('Item')
ylabel('Posterior RT')
title(sprintf('Subject %s', subjects{fitsub}))
legend(vh, 'Control', 'Correlated', 'Filtering')

%
itemNames = {'Near', 'Middle', 'Far'};
fprintf('          \t %14s\t %14s\t %14s\n', itemNames{1}, itemNames{2}, itemNames{3}); 
fprintf(sprintf('Subject %s\n', subjects{fitsub}))
taskNames = {'Control', 'Correlated', 'Filtering'};
for i = 1:3
        fprintf('%10s\t %3.0f (%3.0f, %3.0f)\t %3.0f (%3.0f, %3.0f)\t %3.0f (%3.0f, %3.0f)\n',...
            taskNames{i}, means(i,1), cis(1,1,i), cis(2,1,i),...
                          means(i,2), cis(1,2,i), cis(2,2,i),...
                          means(i,3), cis(1,3,i), cis(2,3,i));

end
% save(sprintf('SEP_%s.mat', subjects{fitsub}))
end
% % cis = prctile(muSamp, [2.5, 97.5]);
% cis = prctile(muSamp, [.5, 99.5]);
% for i = 1:3
% %     text(i+.1, cis(2,i), ['95% HDIs = ', sprintf('\n[%3.2f, %3.2f]', cis(1,i), cis(2,i))])
%     text(i+.1, cis(2,i), ['99% HDIs = ', sprintf('\n[%3.2f, %3.2f]', cis(1,i), cis(2,i))])
% end
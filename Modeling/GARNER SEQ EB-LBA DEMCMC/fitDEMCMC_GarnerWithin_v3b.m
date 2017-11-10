% DE-MCMC Code for modified Garner task data (Little, Wang & Nosofsky, 2016)
%
% Fit all subjects from both conditions (bri/sat) using sequential EB
% model with full parameterization (v4)
%
% theta ~ phi
% data ~ theta
clear variables
clc
close all
format short g

% Set up parallel processing here if necesary (it seems to slow this code down)
parpool('local', 8)
fn = 'GarnerWithin_Bri_Sat_meancentered.mat';

%% Load an organize data
[model, data, n, subnums] = loadData;
% [model, data, n, subnums] = loadDataSep; 
%[model, data, n, subnums] = loadDataBoxcar; 
% [model, data, n, subnums] = loadDataModelRecovery;

%% Parameters for inference
[parms, hyperparms, phiprior] = loadParmSettings(n); % Previously [parms, hyperparms, thetaprior, phiprior], but thetaprior can easily constructed from hyperparms structure
% Returns parm structure with initial values and hypers matrix of hyperparameter values
names = fieldnames(parms);
hnames = fieldnames(hyperparms);

%% DE-MCMC parameters
n.subLevParms = numel(names);
n.hyperParms = numel(hnames);
n.parms  = n.subLevParms + n.hyperParms;
n.chains = 3 * n.subLevParms; % Rule of thumb: 3 x number of subject level parameters or more

n.initburn    = 2500; % 500;    % Number of initial burn in samples
n.migration   = 750; %750;    % Number of iterations that use migration
n.postmigburn = 500; %500;  % Number of post-migration burn in samples without migration
n.mcsamples   = 2000; %2000;   % Number of post-burn in samples
n.burnin      = n.initburn + n.migration + n.postmigburn; % round(n.mc/3); % Guess for burn-in level.
n.migrationStep = 20; %20; % How often (in trials) to run migration (i.e., every x trials)
n.mc          = n.burnin + n.mcsamples; % Total number of samples
use.accept    = 0;

%%
beta     = .001; % Small noise value to add to suggestions
weight = -Inf * ones(n.chains, n.mc, n.subjects); % Initialize likelihood weight vectors for each chain

% Pre-allocate chains
for i = 1:n.subLevParms
    theta.(names{i}) = cell(1,n.subjects);
    theta.(names{i})(:) = {nan(n.chains, n.mc)};
end

%% Generate some fairly widely-spread random start points
% Set up hyperprior chains (phi)
makehpdist = @(hname)(makePriorDistribution(hname, phiprior(strcmp(hnames, hname), :)));  % Cell function for generating distribution objects using settings in phiprior
phipd = cellfun(makehpdist, hnames, 'UniformOutput', false);                              % Apply cell function to make distributions for each hyperparameter
generatePhiValue = @(pdf_object)(random(pdf_object, n.chains, 1));                        % Cell function for generating random values from a distribution object
phiSamples = cellfun(generatePhiValue, phipd, 'UniformOutput', false);                    % Apply cell function to generate random values from each phi distribution
appendNaN = @(x)([x nan(size(x,1), n.mc-1)]);                                             % Cell function to preallocate remaining chains
phi = cell2struct(cellfun(appendNaN, phiSamples, 'UniformOutput', false), hnames);        % Apply cell function to preallocate all chains. First iteration is a sample, rest are nan

% Set up prior distribution parameters (theta)
thetaprior = reshape(cell2mat(struct2cell(hyperparms))', 2, n.subLevParms)';              % Hyperprior values over theta
makedist = @(name)(makePriorDistribution(name, thetaprior(strcmp(names, name), :)));      % Cell function for generating distribution objects
pd = cellfun(makedist, names, 'UniformOutput', false);                                    % Apply cell function to make prior distribution objects for each parameter in theta
generateValue = @(pdf_object)(random(pdf_object));                                        % Cell function for generating random values from a distribution object

% Get individual subject data
sdata = cell(n.subjects, 1);
for k = 1:n.subjects
    % Get subject data for computing likelihood
    sdata{k} = struct('data', data.data(data.data(:,strcmp(data.cols, 'sub')) == subnums(k), :),...
        'cols', {data.cols}, 'stimloc', {data.stimloc(k)}, 'model', data.model);
end

tic
for i = 1:n.chains
    while any(squeeze(weight(i,1,:)) == -Inf) % Ensure valid parms
        for k = 1:n.subjects
            subparms = cellfun(generateValue, pd);          % Sample parameters from the prior for each subject
            parms = cell2struct(num2cell(subparms), names); % Allocate to parm structure

            % Evaluate likelihood for subject k on chain i
            weight(i,1,k) = logDensLike(parms, sdata{k}); % likelihood
            
            if weight(i,1,k) ~= -Inf; % If likelihood is ok, update parms in theta structure
                % Store all parameters for subject k on chain i
                for j = 1:n.subLevParms
                    theta.(names{j}){k}(i,1) = subparms(j); 
                end
            end
        end
    end
end

%% Timing diagnostics
t = nan(n.mc, 1); t(1) = toc; % Store time information for diagnostics

%% Load saved samples if exists
field = cell2struct(hnames, hnames);
if exist(fn, 'file') == 2 % If the file exists
    load(fn, 'i', 'use', 'theta', 'phi', 'weight', 't'); % Load these variables
    starti = i;                                          % Start at i (and sample to n.mc)
else
    starti = 2;                                          % Otherwise, start at 2
end

%% Run the sampler
for i = starti:n.mc
    fprintf('Iteration %d, Time = %3.2f secs\n', i, t(i-1)); tic; % Show the iteration number and restart the timer
    
    % Get previous iterating samples
    use.theta = getChain(theta, i-1, [], 'all');  % Get previous subject level samples
    use.phi = getChain(phi, i-1, [], 'all');      % Get previous hyper parameter samples
    use.like = squeeze(weight(:,i-1,:));          % Get the previous weight
    
    %% Block update hyperparameters
    % If the iteration is migration iteration
    if i >= n.initburn && i < (n.burnin - n.postmigburn) && mod(i, n.migrationStep) == 0
        % Migrate the hyperparameters in pairs
        cnt = 1; % Index of current pair
        for k = 1:n.hyperParms/2 % Cycle through all pairs
            use = hypermigration(use, cnt, phiprior, n, beta); % Hyper migration step
            cnt = cnt + 2;       % Increment pair counter
        end
    else
        cnt = 1;
        for k = 1:n.hyperParms/2
            use = hypercrossover_v2(cnt, use, phiprior, n,  beta); % Sample new values for each chain
            cnt = cnt + 2;
        end
    end
    
    % Do migration or do cross over
    %     fprintf('Updating Subject Level Parms\n')
    if i >= n.initburn && i < (n.burnin - n.postmigburn) && mod(i, n.migrationStep) == 0
        use = migration(use, data, n, beta);         % Migration step
        for k = 1:n.subjects
            for j = 1:n.chains
                weight(j,i,k) = use.like(j,k);
            end
        end
    else
        use = crossover_v3(use, sdata, n, beta);
        weight(:,i,:) = reshape(use.like, [n.chains, 1, n.subjects]);
    end
    
    % Save samples
    theta = getChain(use.theta, i, [], 'update', theta); % Keep new samples
    phi = getChain(use.phi, i, [], 'update', phi);     % Keep new hyper samples
    
    t(i) = toc;
    if mod(i, 200) == 0
        save(fn)
    end
end
totalTime = sum(t);
save(fn)
delete(gcp('nocreate'))
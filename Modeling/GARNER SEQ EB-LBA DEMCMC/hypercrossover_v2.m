function [use] = hypercrossover_v2(kdx, use, phiprior, n, beta)
% Do crossover steps for hyperparameters only
% v4 - do hyperparameter crossover in pairs
%
% idx   = current chain
% kdx   = current hyperparameter
% use   = structure containing sampled parameters
% hyper = structure containing hyper prior values
% n     = structure containing number of chains, mc samples, and parms
% beta  = small value for uniform random noise

% Get parameters from chain idx
names = fieldnames(use.theta);
% theta = getChain(use.theta, idx, [], 'one'); % Chain idx only
% theta = rmfield(theta, names(~ismember(names, names{ceil(kdx./2)}))); % Remove all fields but current parameter 
theta = cell2mat(use.theta.(names{ceil(kdx./2)})); % Get all chains of current parm

getTwoFields = @(x, two)([x.(two{1}), x.(two{2})]);
hnames = fieldnames(use.phi);
% phi   = getChain(use.phi, idx, [], 'one');
% phi = rmfield(phi, hnames(~ismember(hnames, hnames(kdx:kdx+1))));
oldphi = getTwoFields(use.phi, hnames(kdx:kdx+1)); % Get all chains
phi = oldphi;

% Update subject parameter hyperparms based on current chain value
% Relies on structure order being subject level parms followed by the two
% hyper parms for each subject level parm (if different would need to
% update hyperidx)
% thetaprior = cell2mat(struct2cell(phi))';

% Current sample weights
% use.weight(idx) = use.like(idx) + logDensPrior(theta, thetaprior) + logDensPrior(phi, phiprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior x hyperprior)
% Weight for the current sample is the target probability of that sample (i.e., likelihood x prior x hyperprior)
oldweight = logDensPriorAll(theta, phi, names{ceil(kdx./2)}) +...
            logDensPriorPhi(phi, phiprior(kdx:kdx+1,:), hnames(kdx:kdx+1)); % Weight for each chain

% DE-Proposals
gamma = 2.38/sqrt(2 * length(hnames)); % Tuning parameters
chainIdx = 1:n.chains;
selectOtherChains = @(x, i)(datasample(x(x ~= i), 2, 'Replace', false)); % Function for selecting two chains other than the one indexed by i
index = cell2mat(arrayfun(@(x)selectOtherChains(chainIdx, x), chainIdx, 'UniformOutput', false)'); % Array function for applying selectOtherChains to all indexes
% chainIdx(chainIdx == idx) = [];
% index = datasample(chainIdx, 2, 'Replace', false); % Select two other chains

% hnames = hnames(kdx:kdx+1);
% for k = 1:2 %kdx:kdx + 1;
%     phi.(hnames{k}) = use.phi.(hnames{k})(idx,:) +... % Current chain
%         gamma * (use.phi.(hnames{k})(index(1)) - use.phi.(hnames{k})(index(2))) + ... % gamma-scaled difference between two other chains
%         unifrnd(-beta, beta);
% end

% Update phi
phi = phi + phi(index(:,1), :) - phi(index(:,2), :) + unifrnd(-beta, beta, size(phi));
newweight = logDensPriorAll(theta, phi, names{ceil(kdx./2)}) +...
            logDensPriorPhi(phi, phiprior(kdx:kdx+1,:), hnames(kdx:kdx+1)); % Weight for each chain
newweight(any(phi < 0, 2) | any(isnan(phi), 2)) = -Inf;

% check = struct2cell(phi);
% check = [check{:}];

% if checkParms(phi, check, hnames)
%     newweight = -Inf;
% else
%     thetaprior = cell2mat(struct2cell(phi))';
%     thetaprior = reshape(cell2mat(struct2cell(phi))', 2, n.subLevParms)';
%     newweight = use.like(idx) + logDensPrior(theta, thetaprior) + logDensPrior(phi, phiprior); % New weight from like (don't need to recompute as subject level parms haven't changed) and new prior
%     newweight = logDensPrior(theta, thetaprior) + logDensPrior(phi,  phiprior(kdx:kdx+1,:)); % New weight from like (don't need to recompute as subject level parms haven't changed) and new prior
% end

% if isnan(newweight);
%     newweight = -Inf; % Replace any nan weights with -Inf
% end

selectedSet = rand(size(phi,1), 1) < exp(newweight-oldweight);
oldphi(selectedSet,:) = phi(selectedSet,:);
use.phi.(hnames{kdx}) = oldphi(:,1);
use.phi.(hnames{kdx+1}) = oldphi(:,2);
% if rand < (exp(newweight - oldweight)) % rand < new/old, then update theta and like
%     use.phi.(hnames{1})(idx,:) = phi.(hnames{1});
%     use.phi.(hnames{2})(idx,:) = phi.(hnames{2});
% end
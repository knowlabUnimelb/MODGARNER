function use = crossover(idx, kdx, use, data, n, beta, updateChain)
% Crossover subject level parms conditioned on hyperparameters one subject
% at a time
%
% idx   = current chain
% kdx   = current subject
% use   = structure containing sampled parameters
% data  = structure containing rt and acc data
% n     = structure containing number of chains, mc samples, and parms
% beta  = small value for uniform random noise

% Get parameters from chain idx

temp = getChain(use.theta, idx, [], 'one');
names = fieldnames(temp);
subparms = cell2mat(struct2cell(temp));
theta = cell2struct(mat2cell(subparms(:,kdx), ones(size(subparms,1),1), 1), names);

phi   = getChain(use.phi, idx, [], 'one');

% Update subject parameter hyperparms based on current chain value
% Relies on structure order being subject level parms followed by the two
% hyper parms for each subject level parm (if different would need to
% update hyperidx)
thetaprior = reshape(cell2mat(struct2cell(phi))', 2, n.subLevParms)';

% Current sample weights
oldweight = use.like(idx, kdx) + logDensPrior(theta, thetaprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)

% DE-Proposals
gamma = 2.38/sqrt(2 * length(names)); % Tuning parameters
chainIdx = 1:n.chains;
chainIdx(chainIdx == idx) = [];

% Make DE-proposals
index = datasample(chainIdx, 2, 'Replace', false); % Select two other chains
for i = 1:n.subLevParms
    current = use.theta.(names{i}){kdx}(idx,:);
    prop1 = use.theta.(names{i}){kdx}(index(1));
    prop2 = use.theta.(names{i}){kdx}(index(2));
    
    theta.(names{i}) = current + gamma * (prop1 - prop2) + unifrnd(-beta, beta);
end

check = struct2cell(theta);
check = [check{:}];

if checkParms(theta, check, names)
    newweight = -Inf;
else
    like = logDensLike(theta, data);               % New likelihood
    
    phi   = getChain(use.phi, updateChain, [], 'one');
    thetaprior = reshape(cell2mat(struct2cell(phi))', 2, n.subLevParms)';
    newweight = like + logDensPrior(theta, thetaprior);% + logDensPrior(phi, phiprior); % New weight from new like and new prior
end

if isnan(newweight);
    newweight = -Inf; % Replace any nan weights with -Inf
end

if rand < (exp(newweight - oldweight)) % rand < new/old, then update theta and like
    use.accept = use.accept + 1;
    for i = 1:n.subLevParms
        use.theta.(names{i}){kdx}(idx,:) = theta.(names{i});
    end
    use.like(idx,kdx) = like;
end
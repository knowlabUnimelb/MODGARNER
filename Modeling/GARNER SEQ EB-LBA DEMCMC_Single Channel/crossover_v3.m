function use = crossover_v3(use, sdata, n, beta)
% Crossover subject level parms conditioned on hyperparameters 

% use   = structure containing sampled parameters
% sdata  = structure containing rt and acc data
% n     = structure containing number of chains, mc samples, and parms
% beta  = small value for uniform random noise

% Get field names
names = fieldnames(use.theta); 

% Extract theta from use.theta
theta = reshape(cell2mat(cellfun(@(x)(cell2mat(x)), struct2cell(use.theta), 'UniformOutput', false)), [n.chains, n.subLevParms, n.subjects]);
% Extract phi from use.phi
phi   = reshape(cell2mat(struct2cell(use.phi)), n.chains, n.hyperParms); 

% Sampling
gamma = 2.38/sqrt(2 * length(names)); % Tuning parameters
chainIdx = 1:n.chains;
selectOtherChains = @(x, i)(datasample(x(x ~= i), 2, 'Replace', false)); % Function for selecting two chains other than the one indexed by i

% Preallocate
oldweight = nan(size(use.like));
newweight = nan(size(oldweight));
newlike = nan(size(oldweight));

oldtheta = theta;
outTheta = theta;

parfor i = 1:n.subjects
% for i = 1:n.subjects
    % Create large matrices of nchains by nparms for each subject (for theta, phi, and the name of the parameter)
    subTheta = reshape(theta(:,:,i), n.chains * n.subLevParms, 1);
    subPhi = [reshape(phi(:,1:2:end), n.chains * n.subLevParms, 1),... % First hyperparm
              reshape(phi(:,2:2:end), n.chains * n.subLevParms, 1)];   % Second hyperparm
    subParmName = reshape(repmat(names', n.chains, 1), n.subLevParms * n.chains, 1);       
    getIndexedPriorDensity = @(n,x,h,i)getPriorDensity(n{i}, x(i), h(i,:));         % Extract density for a given row
    
    ldp = arrayfun(@(i)getIndexedPriorDensity(subParmName, subTheta, subPhi, i), 1:(n.chains * n.subLevParms))';    % Log of Prior density for each subject's parameters
    oldweight(:,i) = use.like(:,i) + sum(reshape(ldp, n.chains, n.subLevParms), 2); % Weight for each chain and subject

    index = cell2mat(arrayfun(@(x)selectOtherChains(chainIdx, x), chainIdx, 'UniformOutput', false)'); % Array function for applying selectOtherChains to all indexes
    current = theta(:,:,i);
    prop1 = current(index(:,1), :);
    prop2 = current(index(:,2), :);
    current = current + gamma * (prop1 - prop2) + unifrnd(-beta, beta, n.chains, n.subLevParms);

    newSubTheta = reshape(current, n.chains * n.subLevParms, 1);
    newldp = arrayfun(@(i)getIndexedPriorDensity(subParmName, newSubTheta, subPhi, i), 1:(n.chains * n.subLevParms))';    % Log of Prior density for each subject's parameters
    
    likfun = @(j)logDensLike(cell2struct(num2cell(current(j,:))', names), sdata{i});
    newlike(:,i) = cell2mat(arrayfun(@(x)likfun(x), 1:n.chains, 'UniformOutput', false))';
    newweight(:,i) = newlike(:,i) +...
                      + sum(reshape(newldp, n.chains, n.subLevParms), 2);
    
    outTheta(:,:,i) = current;              
end
newweight(squeeze(any(outTheta < 0, 2))) = -Inf; % Keep parms > 0
newweight(squeeze(any(outTheta(:,mstrfind(names, {'w', 'wc', 'b1', 'b2', 'b3', 'b4', 'b5', 'wp'}), :) > 1, 2))) = -Inf; % Keep some parms < 1
newweight(isnan(newweight)) = -Inf;

selectedSet = rand(n.chains, n.subjects) < exp(newweight-oldweight);
use.like(selectedSet) = newlike(selectedSet);
for i = 1:n.subjects
    oldtheta(selectedSet(:,i), :, i) = outTheta(selectedSet(:,i), :, i);
    
    for j = 1:n.subLevParms
        use.theta.(names{j}){i} = oldtheta(:,j,i);
    end
end
function use = crossover_v2(use, sdata, n, beta)
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

% temp = getChain(use.theta, idx, [], 'one');
names = fieldnames(use.theta);
% subparms = cell2mat(struct2cell(temp));
% theta = cell2struct(mat2cell(subparms(:,kdx), ones(size(subparms,1),1), 1), names);
% phi   = getChain(use.phi, idx, [], 'one');

% nchains x nparms x nsubjects matrix
THETA = reshape(cell2mat(cellfun(@(x)(cell2mat(x)), struct2cell(use.theta), 'UniformOutput', false)), [n.chains, n.subLevParms, n.subjects]);
% nchains x nhyperparms
PHI   = reshape(cell2mat(struct2cell(use.phi)), n.chains, n.hyperParms); 

% Update subject parameter hyperparms based on current chain value
% Relies on structure order being subject level parms followed by the two
% hyper parms for each subject level parm (if different would need to
% update hyperidx)
% thetaprior = reshape(cell2mat(struct2cell(phi))', 2, n.subLevParms)';

% Current sample weights
% oldweight = use.like(idx, kdx) + logDensPrior(theta, thetaprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)

gamma = 2.38/sqrt(2 * length(names)); % Tuning parameters
chainIdx = 1:n.chains;
selectOtherChains = @(x, i)(datasample(x(x ~= i), 2, 'Replace', false)); % Function for selecting two chains other than the one indexed by i

oldweight = nan(size(use.like));
newweight = nan(size(oldweight));
newlike = nan(size(oldweight));

index = nan(n.chains, 2, n.subjects);
OLDTHETA = THETA;
for i = 1:n.subjects
    % Create large matrices of nchains by nparms for each subject (for theta, phi, and the name of the parameter)
    subTheta = reshape(THETA(:,:,i), n.chains * n.subLevParms, 1);
    subPhi = [reshape(PHI(:,1:2:end), n.chains * n.subLevParms, 1),... % First hyperparm
              reshape(PHI(:,2:2:end), n.chains * n.subLevParms, 1)];   % Second hyperparm
    subParmName = reshape(repmat(names', n.chains, 1), n.subLevParms * n.chains, 1);       
    getIndexedPriorDensity = @(n,x,h,i)getPriorDensity(n{i}, x(i), h(i,:));         % Extract density for a given row
    
    ldp = arrayfun(@(i)getIndexedPriorDensity(subParmName, subTheta, subPhi, i), 1:(n.chains * n.subLevParms))';    % Log of Prior density for each subject's parameters
    oldweight(:,i) = use.like(:,i) + sum(reshape(ldp, n.chains, n.subLevParms), 2); % Weight for each chain and subject

    index(:,:,i) = cell2mat(arrayfun(@(x)selectOtherChains(chainIdx, x), chainIdx, 'UniformOutput', false)'); % Array function for applying selectOtherChains to all indexes
    current = THETA(:,:,i);
    prop1 = THETA(index(:,1,i), :, i);
    prop2 = THETA(index(:,2,i), :, i);
    THETA(:,:,i)= current + gamma * (prop1 - prop2) + unifrnd(-beta, beta, n.chains, n.subLevParms);

    newSubTheta = reshape(THETA(:,:,i), n.chains * n.subLevParms, 1);
    newldp = arrayfun(@(i)getIndexedPriorDensity(subParmName, newSubTheta, subPhi, i), 1:(n.chains * n.subLevParms))';    % Log of Prior density for each subject's parameters
    
    likfun = @(j)logDensLike(cell2struct(num2cell(THETA(j,:,i))', names), sdata{i});
    newlike(:,i) = cell2mat(arrayfun(@(x)likfun(x), 1:n.chains, 'UniformOutput', false))';
    newweight(:,i) = newlike(:,i) +...
                      + sum(reshape(newldp, n.chains, n.subLevParms), 2);
end
newweight(squeeze(any(THETA < 0, 2))) = -Inf; % Keep parms > 0
newweight(squeeze(any(THETA(:,mstrfind(names, {'w', 'wc', 'b1', 'b2', 'b3', 'b4', 'b5', 'wp'}), :) > 1, 2))) = -Inf; % Keep some parms < 1

selectedSet = rand(n.chains, n.subjects) < exp(newweight-oldweight);
use.like(selectedSet) = newlike(selectedSet);
for i = 1:n.subjects
    OLDTHETA(selectedSet(:,i), :, i) = THETA(selectedSet(:,i), :, i);
    
    for j = 1:n.subLevParms
        use.theta.(names{j}){i} = OLDTHETA(:,j,i);
    end
end



% if checkParms(theta, check, names)
%     newweight = -Inf;
% else
%     like = logDensLike(theta, data);               % New likelihood
%     
%     phi   = getChain(use.phi, updateChain, [], 'one');
%     thetaprior = reshape(cell2mat(struct2cell(phi))', 2, n.subLevParms)';
%     newweight = like + logDensPrior(theta, thetaprior);% + logDensPrior(phi, phiprior); % New weight from new like and new prior
% end

% if isnan(newweight);
%     newweight = -Inf; % Replace any nan weights with -Inf
% end

% % if rand < (exp(newweight - oldweight)) % rand < new/old, then update theta and like
%     use.accept = use.accept + 1;
%     for i = 1:n.subLevParms
%         use.theta.(names{i}){kdx}(idx,:) = theta.(names{i});
%     end
%     use.like(idx,kdx) = like;
% end
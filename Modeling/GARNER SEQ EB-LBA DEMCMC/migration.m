function use = migration(use, data, n, beta)

% Migration set
chainIdx = 1:n.chains;
csetsize = round(n.chains/2);
chainSet = datasample(chainIdx, csetsize, 'Replace', false);
propSet  = circshift(chainSet, [0, 1]);

% Get parameter samples
names = fieldnames(use.theta);
hnames = fieldnames(use.phi);

chain(n.chains).theta = [];
chain(n.chains).phi   = [];
for idx = 1:csetsize
    for i = 1:n.subLevParms;
        for k = 1:n.subjects
            chain(chainSet(idx)).theta.(names{i})(k)     = use.theta.(names{i}){k}(chainSet(idx));                                    % Get parms from current chain
            chain(chainSet(idx)).proptheta.(names{i})(k) = use.theta.(names{i}){k}(propSet(idx)) + unifrnd(-beta, beta);             % Get parms from current chain
        end
    end
    
    for i = 1:n.hyperParms
        chain(chainSet(idx)).phi.(hnames{i})     = use.phi.(hnames{i})(chainSet(idx));                                    % Get parms from current chain
    end
    thetaprior = reshape(cell2mat(struct2cell(chain(chainSet(idx)).phi))', 2, n.subLevParms)';
    
    % Current sample weights
    chain(chainSet(idx)).oldweight = sum(use.like(chainSet(idx),:)) +...
        logDensPrior(chain(chainSet(idx)).theta, thetaprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
    
    % Proposal samples weights
    check = struct2cell(chain(chainSet(idx)).proptheta); check = [check{:}];
    if checkParms(chain(chainSet(idx)).proptheta, check, names) % && checkParms(chain(idx).propphi, checkphi, hnames)
        chain(chainSet(idx)).newweight = -Inf;
    else
        subnums = unique(data.data(:, strcmp(data.cols, 'sub')));
        subparms = cell2mat(struct2cell(chain(chainSet(idx)).proptheta));
        for i = 1:n.subjects
            sdata = struct('data', data.data(data.data(:,strcmp(data.cols, 'sub')) == subnums(i), :),...
                'cols', {data.cols}, 'stimloc', {data.stimloc(i)}, 'model', data.model);
            theta = cell2struct(mat2cell(subparms(:,i), ones(size(subparms,1),1), 1), names);
            
            check = struct2cell(theta);
            check = [check{:}];
            if checkParms(theta, check, names)
                lik(i) = -Inf;
            else
            	lik(i) = logDensLike(theta, sdata);
            end
        end
        
        chain(chainSet(idx)).newlike = lik; % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
        chain(chainSet(idx)).newweight =  sum(chain(chainSet(idx)).newlike) +...
           logDensPrior(chain(chainSet(idx)).proptheta, thetaprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
    end
end

% Migration
for i = 1:csetsize
    if rand < (exp(chain(chainSet(i)).newweight - chain(chainSet(i)).oldweight))
        fprintf('Migrating...\n')
        for j = 1:numel(names)
            for k = 1:n.subjects
                use.theta.(names{j}){k}(chainSet(i)) = chain(chainSet(i)).proptheta.(names{j})(k);
            end
        end
        use.like(chainSet(i), :) = chain(chainSet(i)).newlike;
    end
end
% HYPERMIGRATION
%
% Migrate hyperparms in pairs (Turner, 2013, p. 384)
%
% Propose a jump from one chain’s current state to another chain’s current 
% state. This proposal often includes multiple chain states being swapped 
% in a cyclical fashion, so that Chain 1 moves to Chain 2 and Chain 2 moves 
% to Chain 3 and Chain 3 moves to Chain 1. 
function use = hypermigration(use, kdx, phiprior, n, beta)

% Migration set
chainIdx = 1:n.chains;                                       % Index for each chain
csetsize = round(n.chains/2);                                % Swap half of chains
chainSet = datasample(chainIdx, csetsize, 'Replace', false); % Sample the indexes without replacement
propSet  = circshift(chainSet, [0, 1]);                      % Shift the sampled indexes by 1 to form the proposal indexes

% Get parameter samples
names = fieldnames(use.theta); % Theta parm names
hnames = fieldnames(use.phi);  % Phi parm names

groupCells = @(x)([x{:}]);

subjectTheta = groupCells(use.theta.(names{ceil(kdx/2)})); % Get all chains x all subjects for the theta corresponding to kdx
subjectThetaChains = subjectTheta(chainSet,:);             % Grab only the chains being migrated

phi = cell2mat(struct2cell(use.phi));        % Extract phi from structure
phi = reshape(phi, n.chains, numel(hnames)); % Reshape phi
phi = phi(:,kdx:kdx+1);                      % Keep only selected parms
current = phi(chainSet, :);
proposal = phi(propSet, :) + unifrnd(-beta, beta, numel(propSet), 2);

% Current samples weights for each target chain for given hyperparm
oldweight = logDensPrior_v2(subjectThetaChains, current, names(ceil(kdx/2))) +...
            logDensPrior_v2(current, phiprior, hnames(kdx:kdx+1)); 

% Proposal samples weights for each target chain for given hyperparm        
newweight = logDensPrior_v2(subjectThetaChains, proposal, names(ceil(kdx/2))) +...
            logDensPrior_v2(proposal, phiprior, hnames(kdx:kdx+1));          
newweight(any(proposal < 0, 2)) = -Inf;  % Check that the proposal parms are sensible

updateset = rand(csetsize, 1) < exp(newweight - oldweight); % Flag chains to update
use.phi.(hnames{kdx})(chainSet(updateset)) = proposal(updateset, 1);
use.phi.(hnames{kdx+1})(chainSet(updateset)) = proposal(updateset, 2);
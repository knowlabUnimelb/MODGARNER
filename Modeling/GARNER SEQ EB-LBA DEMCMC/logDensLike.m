function [out, lik] = logDensLike(parms, data)
% Organize parameters based on the model
% Sum these loglikelihoods to compute the overall log likelihood

%% Set up the model 
fmodel = data.model;

%% Set up subject level parms
iparms = cell2mat(struct2cell(parms))';

%% Loop
subs = unique(data.data(:, strcmp(data.cols, 'sub')));
lik = nan(numel(subs),1);
% parfor i = 1:numel(subs)
for i = 1:numel(subs)
    % Get likelihood for each subject
    lik(i) = fmodel(iparms(i,:), data.stimloc{i},...
        data.data(data.data(:,strcmp(data.cols, 'sub')) == subs(i), :),...
        data.cols);
end
out = sum(lik);
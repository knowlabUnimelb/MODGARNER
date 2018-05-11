function out = getChain(in, idx, n, varargin)

optargs = {'one', [], []};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[whichfun, out, beta] = optargs{:}; % Place optional args in memorable variable names

%%
switch whichfun
    case 'one'
        % Get parameters from chain idx
        names = fieldnames(in);
        switch class(in.(names{1}))
            case 'cell' % Subject level parms
                out = structfun(@(x) getOneTheta(x, idx), in, 'UniformOutput', false);
            case 'double'
                out = structfun(@(x) getOnePhi(x, idx), in, 'UniformOutput', false);
        end
    case 'all'
        % Get all chains for given iteration
        names = fieldnames(in);
        switch class(in.(names{1}))
            case 'cell' % Subject level parms
                for j = 1:numel(names)
                    for k = 1:size(in.(names{j}), 2)
                        out.(names{j}){k} = in.(names{j}){k}(:,idx); % Get the indicated sample
                    end
                end
            case 'double' % Hyper parms
                for j = 1:numel(names)
                    out.(names{j}) = in.(names{j})(:,idx); % Get the indicated sample
                end
        end
    case 'update'
        names = fieldnames(in);
        switch class(in.(names{1}))
            case 'cell' % Subject level parms
                for j = 1:numel(names)
                    for k = 1:size(in.(names{j}), 2)
                        out.(names{j}){k}(:,idx) = in.(names{j}){k};   % Assign samples to sampled parameter chains
                    end
                end
            case 'double' % Hyper parms
                for j = 1:numel(names)
                    out.(names{j})(:,idx) = in.(names{j});   % Assign samples to sampled parameter chains
                end
        end
    case 'prop'
        % Get parameters from chain idx
        names = fieldnames(in);  % Get parameter names
        for i = 1:numel(names);         % Cycle through parms
            if i <= n.subLevParms       % If subject level parms, then retrive from each cell
                for k = 1:n.subjects
                    out.(names{i})(1,k) = in.(names{i}){k}(idx) + unifrnd(-beta, beta);
                end
            else                        % Otherwise, just get the parm
                out.(names{i}) = in.(names{i})(idx) + unifrnd(-beta, beta);
            end
        end
end


function out = getOneTheta(in, idx)

for k = 1:size(in, 2)
    out(1,k) = in{k}(idx);
end

function out = getOnePhi(in, idx)
out = in(idx);

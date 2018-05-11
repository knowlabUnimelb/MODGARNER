function theta = logitParms(theta, direction)

switch direction
    case 'unbound'
        names = fieldnames(theta);
        for i = 1:numel(names)
            if ismember(names{i}, {'w', 'wc', 'b1', 'b2', 'b3', 'b4', 'b5'})
                theta.(names{i}) = logit(theta.(names{i}));
            end
        end    
    case 'bound'
        names = fieldnames(theta);
        for i = 1:numel(names)
            if ismember(names{i}, {'w', 'wc', 'b1', 'b2', 'b3', 'b4', 'b5'})
                theta.(names{i}) = logit(theta.(names{i}), 'inverse');
            end
        end
end
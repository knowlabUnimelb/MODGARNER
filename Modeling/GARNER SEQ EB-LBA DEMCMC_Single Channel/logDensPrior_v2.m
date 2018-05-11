function out = logDensPrior_v2(x, hypers, name)
% v2 - specify name as input

out = zeros(size(x,1), 1);
for i = 1:numel(name)
    for j = 1:size(x,2)
        pdf = getPriorDensity(name{i}, x(:,j), hypers); % Get density in log
        out = out + sum(pdf, 2);
    end
end
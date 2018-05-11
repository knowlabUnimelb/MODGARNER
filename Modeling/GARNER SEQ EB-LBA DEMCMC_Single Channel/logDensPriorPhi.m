function out = logDensPriorPhi(x, hypers, name)
% Compute lnL over all subjects

out = zeros(size(x,1), 1);
for i = 1:numel(name)
    pdf = getPhiPriorDensity(name{i}, x(:,i), hypers(i,:)); % Get density in log
    out = out + sum(pdf, 2);
end

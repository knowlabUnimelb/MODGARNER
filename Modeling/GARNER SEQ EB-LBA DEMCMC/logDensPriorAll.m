function out = logDensPriorAll(x, hypers, name)
% Compute lnL over all subjects

out = zeros(size(x,1), 1);
for i = 1:size(x,2)
    pdf = getPhiPriorDensity(name, x(:,i), hypers); % Get density in log
    out = out + sum(pdf, 2);
end

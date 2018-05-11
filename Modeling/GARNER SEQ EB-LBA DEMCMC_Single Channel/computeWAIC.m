function WAIC = computeWAIC(data, theta, model, n)

names = fieldnames(theta);
psamples = round(.1 * n.mcsamples);

for midx = 1:psamples
    for iidx = 1:numel(names);
        use.theta.(names{iidx}) = theta.(names{iidx})(datasample(1:n.chains, 1), datasample(n.burnin+1:n.mcsamples, 1));             % Get parms from current chain
    end
    
    [~, lnL] = logDensLikeLR(use.theta, data, model);
    lnl(:,midx) = cell2mat(lnL');
    if mod(midx, 10) == 0
        disp(midx)
    end
end
WAIC = 2 * sum(log(sum(exp(lnl), 2)./psamples) - sum(lnl,2)./psamples);


    
    
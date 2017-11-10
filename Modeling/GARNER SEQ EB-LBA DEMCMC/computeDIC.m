function dic = computeDIC(model, data, theta, weight, nburnin, version)

dev = -2 * weight(:); 
dev_hat = mean(dev); 

avgparms = []; 
names = fieldnames(theta);
for i = 1:numel(names)
        temp = theta.(names{i});
        temp = temp(:,nburnin:end);
        eval(sprintf('%s = mean(temp(:));', names{i}));
        
        eval(sprintf('avgparms = setfield(avgparms, names{i}, %s);', names{i}));
end

% Compute loglikelihood for average sample
if version == 1
avglnL = logDensLikeLR(avgparms, data, model);
dev_avg = -2 * avglnL;
pd = dev_hat - dev_avg;
else
pd = var(dev)/2;
end
dic = dev_hat + pd;
function [dic, dev_avg, pd] = computeDICgroup(data, theta, weight, version, idx)

%% Gelman method
subnums = unique(data.data(:,strcmp(data.cols, 'sub')));
nsubs = numel(subnums);

sdata = cell(nsubs, 1);
for k = 1:nsubs
    % Get subject data for computing likelihood
    sdata{k} = struct('data', data.data(data.data(:,strcmp(data.cols, 'sub')) == subnums(k), :),...
        'cols', {data.cols}, 'stimloc', {data.stimloc(k)}, 'model', data.model);
end


avgparms = [];
names = fieldnames(theta);
for sidx = 1:nsubs
    for i = 1:numel(names)
        %     eval(sprintf('avgparms = setfield(avgparms, names{i}, cell(1,nsubs));'));
        
        tempt = theta.(names{i});
        temp = tempt{sidx}(:,idx);
        eval(sprintf('%s = mean(temp(:));', names{i}));
%         eval(sprintf('%s = mode(temp(:));', names{i}));
        eval(sprintf('avgparms.(names{i}) = %s;', names{i}));
        rndparms.(names{i}) = datasample(temp(:), 1);
    end
    sap(sidx) = avgparms;
    srp(sidx) = rndparms;
    dev_avg(sidx) = logDensLike(avgparms, sdata{sidx});
    dev_rnd(sidx) = logDensLike(rndparms, sdata{sidx});
    tempw = weight(:,:,sidx);
    subweightmean(:,sidx) = tempw(:);
end

if version == 1
   pd = 2 * (dev_avg - mean(subweightmean)); 
   save temp1
else
   pd = 2 * var(subweightmean);
   save temp2
end
dic = -2 * sum(dev_avg - pd);
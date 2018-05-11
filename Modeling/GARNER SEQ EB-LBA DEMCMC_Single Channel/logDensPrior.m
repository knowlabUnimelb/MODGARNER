function out = logDensPrior(x, hypers)

names = fieldnames(x);
out = 0; 
for j = 1:numel(names)
    xval = x.(names{j});
    pdf = getPriorDensity(names{j}, xval, hypers(j,:)); % Get density in log
    out = out + sum(pdf);
end
function out = logDensPriorCrossover(x, hypers, name)
% x is a vector of theta
% hypers is a numel(x) x 2 vector of phi
% names is  vector of names



out = 0; 
for j = 1:numel(names)
    xval = x.(names{j});
    pdf = getPriorDensity(names{j}, xval(j), hypers(j,:)); % Get density in log
    out = out + sum(pdf);
end
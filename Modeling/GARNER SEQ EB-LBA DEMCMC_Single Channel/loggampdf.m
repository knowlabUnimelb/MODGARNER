function y = loggampdf(x,a,b)

x = x(:);

if (a < 0 || b < 0)
    y = -Inf * ones(1, numel(x));
else
    y = nan(numel(x), 1);
    y(x < 0) = -Inf; 
    
    a1 = zeros(numel(x), 1);
    a1(x > 0) = (a - 1) * log(x(x > 0));
    a1(x == 0) = 0;
    
    y(x >= 0) = -gammaln(a) + a .* log(b) + a1(x >= 0) - (x(x >= 0) .* b);
end


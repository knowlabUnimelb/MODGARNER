function y = standardNormalCdf(x)
% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
y = (0.5 * erfc(-x ./ sqrt(2)));
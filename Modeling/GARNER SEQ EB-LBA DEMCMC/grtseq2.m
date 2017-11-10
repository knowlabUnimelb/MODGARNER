function p = test_grtseq(stimulusCoordinates, prev, item, a, b, c, cmod, s, smod, igp)
% GRT Apply decision boundary theory to compute the probability that each
% region of a bivariate normal distribution falls within each category 
%
% The decision bounds are in the form aX + bY + c = 0
%
% x - value of stimulus on dimension x
% y - value of stimulus on dimesnion y
% a - coefficient on x (set to 0 if orthogonal)
% b - coefficient on y (set to 1)
% c - intercept
% 
% In this model, the idea is that:
%   1) If the previous item is the same, then the perceptual variance is
%   multiplied by 1/smod; otherwise, the perceptual variance is increased
%   by multiplying by smod
%   2) If the previous category is category A, then the boundary shifts
%   either toward or away from category A and away from or toward category
%   B, respectively

% dist = a(1) .* x(1) + b .* y(1) + c(1);
% dist = dist./sqrt(a(1).^2 + b.^2);
% z = -dist./s(1); 
% 

ix = stimulusCoordinates(item, 1); 
iy = stimulusCoordinates(item, 2);

px = stimulusCoordinates(prev, 1); 
py = stimulusCoordinates(prev, 2);

cat = ones(numel(item), 1);
cat(item >= 7) = 2;

pcat = ones(numel(prev), 1);
pcat(prev >= 7) = 2;

p = nan(numel(ix), 1);
for i = 1:numel(ix)
    if igp(i) == 1
        sadj = s(i);
    elseif igp(i) ~= 1 && ix(i) == px(i) && iy(i) == py(i) % Previous item is the same
        sadj = s(i) * (1/smod); % Current item variance reduces
    else
        sadj = s(i) * 1; % Current item variance increases
    end
    
    if igp(i) == 1
        cadj = c(i);
    elseif igp(i) ~= 1 && pcat(i) == 1
        cadj = c(i) + cmod; % If cmod is pos, move bound away; if cmod is neg, move bound toward
    else
        cadj = c(i) + -cmod; % If cmod is pos, move bound toward; if cmod is neg, move bound away
    end
    dist = a(i) .* ix(i) + b .* iy(i) + cadj;
    dist = dist./sqrt(a(i).^2 + b.^2);
    z = -dist./sadj; 
    p(i) = normcdf(z); % p = pzscor(z)
end
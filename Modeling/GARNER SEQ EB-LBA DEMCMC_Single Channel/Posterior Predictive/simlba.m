function [rt, resp] = simlba(n, parms)
% Generate random samples from the serial LBA
% parms should include v (1 x 2 A/B accumulators)

names = fieldnames(parms);
for i = 1:numel(names)
    if numel(parms.(names{i})) > 1
        dstr = repmat('%d,', 1, numel(parms.(names{i})));
        eval(sprintf(['%s(:,1) = [' dstr(1:end-1) ']'';'], names{i}, parms.(names{i})));
    else
        eval(sprintf('%s(:,1) = [%d];', names{i}, parms.(names{i})));
    end
end
v = [vc, 1 - vc];
b = A + bMa;

%% Simulate model
nchoices=size(v,2);

k1 = A * rand(n, nchoices); % Random starting points

% Sample drift rates
% Since both drift rates are sampled independently, this can result in both drifts being negative. 
d1 = ones(n,1) * v(1,:) + randn(n, nchoices) * s;

ttfa1            = (b(1)-k1(:,1))./d1(:,1);          % Time is distance divided by rate
ttfa1(ttfa1 < 0)  = Inf;                 % Set any negative times equal to nan
ttfb1            = (b(2)-k1(:,2))./d1(:,2);          % Time is distance divided by rate
ttfb1(ttfb1 < 0)  = Inf;                 % Set any negative times equal to nan
ttf1 = [ttfa1, ttfb1];
[rt, resp] = nanmin(ttf1, [], 2); % Find the minimum time for an accumulator to finish

rt = rt + t0;            % Add t0
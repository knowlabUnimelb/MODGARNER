function varargout = seqEBLBA(parms, stimloc, data, cols, varargin)
% same as seqSTMLBAv4-Garner model for GARNERWITHIN project modified for DE-MCMC

optargs = {false};
newVals = cellfun(@(x) ~isempty(x), varargin);% skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[simulate] = optargs{:}; % Place optional args in memorable variable names

R_METRIC = 2;
SHEPARD_EXPONENT = 2;

%% Set up relevant data 
item = data(:,strcmp(cols, 'item'));
task = data(:,strcmp(cols, 'task'));
cat =  data(:,strcmp(cols, 'cat'));
opp =  data(:,strcmp(cols, 'opp'));
acc =  data(:,strcmp(cols, 'acc'));
rt  =  data(:,strcmp(cols, 'rts'));
exc =  data(:,strcmp(cols, 'exc'));
igp =  data(:,strcmp(cols, 'igp'));

nD   = size(rt, 1);
x = stimloc(item, 1); 
y = stimloc(item, 2);
exemplars = [stimloc, [1 1 1 1 1 1 2 2 2 2 2 2]'];

nchoices = 2;

if size(parms,1) ~= 1; parms = parms'; end

%% seqGCM parameters
nGCMparms = 11;
w1 = parms([1 1 2 2 1]);    % Attention to dimension X
w2 = 1 - w1;
w = [w1; w2];
c = (parms([3 3 3 3 4])); % Specificity
r = repmat(R_METRIC, 1, 5); % r-Metric
b1 = parms(5:9);
b2 = 1 - b1;
br = [b1; b2];
wp = parms(10);                  % Probability of using previous stimulus
bh = parms(11);
shep = SHEPARD_EXPONENT;

% p = stmgcm(exemplars, item, task, w, c, r, br, wp, bh, igp);
p = stmgcm2(exemplars, item, task, w, c, r, br, wp, bh, igp, shep);

vcorr = p; %(itemIdx);

%% LBA parameters
startPoint = (parms(nGCMparms+1)); % Range of uniform start point distribution
threshold  = (parms(nGCMparms+2:nGCMparms+3)) + startPoint; % Accumulator thresholds for each category
vs         = (parms(nGCMparms+4)); % Between-trial drift rate variability
t0         = (parms(nGCMparms+5)); % Nondecision time

%% Set up LBA parameters
vcorr(cat == 2) = 1 - vcorr(cat == 2); % If the correct category = 2, then the drift rate for the correct accumulator is 1 - v
verr = 1 - vcorr;                      % Drift rate for the other category
v = [vcorr, verr]; % Drift rates

bcorr = threshold(cat)';     % Threshold for correct 
berr  = threshold(opp)';     % Threshold for error
b = [bcorr, berr];

%% Use LBA to get the fit of the drift rate
if ~simulate
    vIdx = acc;           % The appropriate drift rate to use is determined by the accuracy
    vIdx(vIdx == 0) = 2;
    vEidx = acc + 1;
    
    idx = (1:nD)';
    idx(vIdx == 2) = idx(vIdx == 2) + nD;
    
    jdx = (1:nD)';
    jdx(vEidx == 2) = jdx(vEidx == 2) + nD;
    
    t = rt - t0;
    
    fi = lbapdfvec(t, v(idx), b(idx), startPoint, vs, 7);
    Fj = lbacdfvec(t, v(jdx), b(jdx), startPoint, vs, 7);
    
    fit = nansum(log(fi(exc == 0)) + log(1 - Fj(exc == 0)));
    
    if ~isreal(fit)
        fit = -realmax;
    end
    varargout{1} = fit;
else
    k1 = startPoint * rand(nD, nchoices); % Random starting points
    
    % Sample drift rates
    % Since both drift rates are sampled independently, this can result in both drifts being negative.
    d1 = v + randn(nD, nchoices) * vs;
    
    % Convert to time
    ttf = (b - k1)./d1;          % Time is distance divided by rate
    ttf(ttf < 0)  = Inf;  % Set any negative times equal to nan
    [rt, acc] = nanmin(ttf, [], 2);
    rt = rt + t0;            % Add t0
    acc(acc == 2) = 0; % Recode as correct/incorrect
    
    varargout{1} = rt;
    varargout{2} = acc;
end

%%
function p = stmgcm2(exemplars, item, task, w, c, r, b, wp, bhigher, igp, shep)
% GCM Apply the generalized context model to compute the probability of
% retrieving exemplars from each category
% 
% Vectorized version 15-Aug-2016 DL
x   = exemplars(:,1);
y   = exemplars(:,2);
cat = exemplars(:,3);

taskItems{1} = 1:2:12;
taskItems{2} = 2:2:12;
taskItems{3} = [1 3 5 8 10 12];
taskItems{4} = [2 4 6 7 9 11];
taskItems{5} = 1:12;

%% Compute long term component
sumSA = nan(12,5);
sumSB = nan(12,5);
ltmp = nan(12,5);
for i = 1:5
    d = ((w(1,i) * (abs(repmat(x,1,numel(x)) - repmat(x', numel(x), 1)).^r(i))) +...
         (w(2,i) * (abs(repmat(y,1,numel(y)) - repmat(y', numel(y), 1)).^r(i)))).^(1/r(i));
    
    s = exp(-c(i) * d.^shep);
    
    ts = s(:,taskItems{i});
    tc = cat(taskItems{i});
    
    sumSA(:,i) = (b(1,i) * sum(ts(:,tc == 1), 2));
    sumSB(:,i) = (b(2,i) * sum(ts(:,tc == 2), 2));
    
    ltmp(:,i) = sumSA(:,i)./(sumSA(:,i) + sumSB(:,i));
end

%% Compute short term component
% Remove first item and treat that item separately
pitem = item(1:end-1);
iitem = item(2:end); % Current item
itask = task(2:end); % Current task

px = exemplars(item(1:end-1),1); % Previous item, dimension 1
ix = exemplars(item(2:end),  1); % Current item, dimension 1 

py = exemplars(item(1:end-1),2); % Previous item, dimension 2 
iy = exemplars(item(2:end),  2); % Current item, dimension 2 

pcat = exemplars(item(1:end-1),3); % Previous category
icat = exemplars(item(2:end),3);   % Current category

boostA = zeros(numel(item)-1,1);  % Preallocate
boostB = zeros(numel(item)-1, 1); % Preallocate

biasA = ones(numel(item)-1, 1); % Preallocate
biasB = ones(numel(item)-1, 1); % Preallocate

simd = (w(1, itask)' .* abs(ix - px).^(r(itask)') + w(2, itask)' .* abs(iy - py).^(r(itask)')).^(1./(r(itask)'));
simprev = exp(-(c(itask)') .* simd.^shep);

boostA(pcat == 1) = wp * simprev(pcat == 1);
boostB(pcat == 2) = wp * simprev(pcat == 2);

biasB(pcat == 1 & iitem ~= pitem) = bhigher;
biasA(pcat == 2 & iitem ~= pitem) = bhigher;

newSA = (biasA .* (sumSA(sub2ind(size(sumSA), iitem, itask)) + boostA));
newSB = (biasB .* (sumSB(sub2ind(size(sumSB), iitem, itask)) + boostB));
p = newSA./(newSA + newSB);

p(igp(2:end) == 1) = ltmp(sub2ind(size(ltmp), iitem(igp(2:end) == 1), itask(igp(2:end) == 1)));
p1 = ltmp(sub2ind(size(ltmp), item(1), task(1)));
p = [p1; p];
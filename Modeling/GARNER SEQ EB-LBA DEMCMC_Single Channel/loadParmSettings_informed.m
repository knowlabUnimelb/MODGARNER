function [parms, hyperparms, phiprior] = loadParmSettings(n)
global SHEPARD_EXPONENT R_METRIC

R_METRIC = 1;
SHEPARD_EXPONENT = 1; % Shepard Exponent (1 = exponential, 2 = Gaussian)
% I want this to be available for editing here

%% EB parms
%% Attention weights
% Attention to irrelevant dimension
w  = .20; % Baseline & Filtering attention
% w ~ Beta(a, b)
w_a = 6; % w_a ~ Normal(1, _5)T(0_5, inf)
w_b = 14; %5; % w_b ~ Normal(1, _5)T(0_5, inf)

w_a_mu    = 6; w_a_sigma = .5;
% w_b_mu    = 1; w_b_sigma = .5;
w_b_mu    = 14; w_b_sigma = .5;

%% 
% Attention to irrelevant dimension
wc = .2; % Correlated condition

% wc ~ Beta(a, b)
wc_a = 5; % w_a ~ Normal(1, _5)T(0_5, inf)
wc_b = 14; %5; % w_b ~ Normal(1, _5)T(0_5, inf)

wc_a_mu    = 5; wc_a_sigma = .5;
% wc_b_mu    = 1; wc_b_sigma = .5;
wc_b_mu    = 14; wc_b_sigma = .5;

%% Sensitivity parameters
c  = 1;  % Baseline & Correlated specificity

% c ~ Normal(mu, sigma)T(0, inf)
c_mu    = 1;
c_sigma = 2;

c_mu_mu = 1; c_mu_sigma = 3;   % Fixed hyper prior values
% c_mu_mu = .5; c_mu_sigma = 1;   % Fixed hyper prior values
c_sigma_a = 1; c_sigma_b = 1;

%%
cf = 1;  % Filtering specificity

cf_mu    = 1;
cf_sigma = 2;

cf_mu_mu = 1; cf_mu_sigma = 3;     % Fixed hyper prior values
% cf_mu_mu = .5; cf_mu_sigma = 1;     % Fixed hyper prior values
cf_sigma_a = 1; cf_sigma_b = 1;

%% Bias parms
b1 = .5; % Baseline 1 - A bias
b2 = .5; % Baseline 2 - A bias
b3 = .5; % Correlated 1 - A bias
b4 = .5; % Correlated 2 - A bias
b5 = .5; % Filtering - A bias

% b1 ~ Beta(a, b)
b1_a = 5; % b1_a ~ Normal(1, _5)T(0_5, inf)
b1_b = 5; % b1_b ~ Normal(1, _5)T(0_5, inf)

b1_a_mu    = 5; b1_a_sigma = .5;
b1_b_mu    = 5; b1_b_sigma = .5;

b2_a = 5; % b2_a ~ Normal(1, _5)T(0_5, inf)
b2_b = 5; % b2_b ~ Normal(1, _5)T(0_5, inf)

b2_a_mu    = 5; b2_a_sigma = .5;
b2_b_mu    = 5; b2_b_sigma = .5;

b3_a = 5; % b3_a ~ Normal(1, _5)T(0_5, inf)
b3_b = 5; % b3_b ~ Normal(1, _5)T(0_5, inf)

b3_a_mu    = 5; b3_a_sigma = .5;
b3_b_mu    = 5; b3_b_sigma = .5;

b4_a = 5; % b4_a ~ Normal(1, _5)T(0_5, inf)
b4_b = 5; % b4_b ~ Normal(1, _5)T(0_5, inf)

b4_a_mu    = 5; b4_a_sigma = .5;
b4_b_mu    = 5; b4_b_sigma = .5;

b5_a = 5; % b5_a ~ Normal(1, _5)T(0_5, inf)
b5_b = 5; % b5_b ~ Normal(1, _5)T(0_5, inf)

b5_a_mu    = 5; b5_a_sigma = .5;
b5_b_mu    = 5; b5_b_sigma = .5;

%% Sequential EB parms
%%
% wp = 1;  % Weight of previous item 
% wp ~ Normal(mu, sigma)T(0, inf)
% 29/8/16 - this is different from the paper where wp is bounded 0, 1
% wp_mu = 1;     % wp_mu ~ Normal(1, _5)T(0, inf)
% wp_sigma = .5; % wp_sigma ~ Gamma(1, 1)
% wp_mu_mu = 1; wp_mu_sigma = .5;
% wp_sigma_a = 1; wp_sigma_b = 1;

wp  = .8; % Baseline & Filtering attention
% wp ~ Beta(a, b)

wp_a = 1; % wp_a ~ Normal(1, _5)T(0_5, inf)
wp_b = 1; % wp_b ~ Normal(1, _5)T(0_5, inf)

wp_a_mu    = 1; wp_a_sigma = .5;
wp_b_mu    = 1; wp_b_sigma = .5;

%%
bh = 1.25;  % Change bias parameter - bhigher(beta)

% bh ~ Normal(mu, sigma)T(0, inf)
bh_mu = 1.25;     % wp_mu ~ Normal(1, _5)T(0, inf)
bh_sigma = 2; % wp_sigma ~ Gamma(1, 1)

bh_mu_mu = 1.25; bh_mu_sigma = 3;
bh_sigma_a = 1; bh_sigma_b = 1;

%% LBA parms
%% Startpoint
A    = .075;  

% A ~ Normal(_12, _5)T(0, inf)
A_mu = .075; 
A_sigma = 2;

A_mu_mu = .075; A_mu_sigma = 3;
A_sigma_a = 1; A_sigma_b = 1;

%% A boundary minus start point
bMa1 = .35; % Threshold minus startpoint

% bMa1 ~ Normal(_2, _5)T(0, inf)
bMa1_mu = .35; 
bMa1_sigma = 2;

bMa1_mu_mu = .35; bMa1_mu_sigma = 3;
bMa1_sigma_a = 1; bMa1_sigma_b = 1;

%% B boundary minus start point
bMa2 = .35; % Threshold minus startpoint

% bMa2 ~ Normal(_2, _5)T(0, inf)
bMa2_mu = .35; 
bMa2_sigma = 2;

bMa2_mu_mu = .35; bMa2_mu_sigma = 3;
bMa2_sigma_a = 1; bMa2_sigma_b = 1;

%% Drift rate variability
vs   = .25;  % Drift variability

% vs ~ Normal(_25, _5)T(0, inf)
vs_mu = .25; 
vs_sigma = 1;

vs_mu_mu = .25; vs_mu_sigma = 3;
vs_sigma_a = 1; vs_sigma_b = 1;

%% Non-decision time
t0   = .1;   % Non-decision time

% t0 ~ Normal(_1, _5)T(0, inf)
t0_mu = .1; 
t0_sigma = 5;

t0_mu_mu = .1; t0_mu_sigma = 10;
t0_sigma_a = 1; t0_sigma_b = 1;

%% Parm structure
% Set up structure which hold the parameters for that model (need
% one section for each model with different parameterization)
parms = struct(...
    'w', w * ones(1,n.subjects),...
    'wc', wc  * ones(1,n.subjects),...
    'c', c  * ones(1,n.subjects),...
    'cf', cf  * ones(1,n.subjects),...
    'b1', b1  * ones(1,n.subjects),...
    'b2', b2  * ones(1,n.subjects),...
    'b3', b3  * ones(1,n.subjects),...
    'b4', b4  * ones(1,n.subjects),...
    'b5', b5  * ones(1,n.subjects),...
    'wp', wp  * ones(1,n.subjects),...
    'bh', bh  * ones(1,n.subjects),...
    'A', A  * ones(1,n.subjects),...
    'bMa1', bMa1  * ones(1,n.subjects),...
    'bMa2', bMa2  * ones(1,n.subjects),...
    'vs', vs  * ones(1,n.subjects),...
    't0', t0  * ones(1,n.subjects));

hyperparms = struct(...
    'w_a', w_a,...
    'w_b', w_b,...
    'wc_a', wc_a,...
    'wc_b', wc_b,...
    'c_mu', c_mu,...
    'c_sigma', c_sigma,...
    'cf_mu', cf_mu,...
    'cf_sigma', cf_sigma,...
    'b1_a', b1_a,...
    'b1_b', b1_b,...
    'b2_a', b2_a,...
    'b2_b', b2_b,...
    'b3_a', b3_a,...
    'b3_b', b3_b,...
    'b4_a', b4_a,...
    'b4_b', b4_b,...
    'b5_a', b5_a,...
    'b5_b', b5_b,...
    'wp_a', wp_a,...
    'wp_b', wp_b,...
    'bh_mu', bh_mu,...
    'bh_sigma', bh_sigma,...
    'A_mu', A_mu,...
    'A_sigma', A_sigma,...
    'bMa1_mu', bMa1_mu,...
    'bMa1_sigma', bMa1_sigma,...
    'bMa2_mu', bMa2_mu,...
    'bMa2_sigma', bMa2_sigma,...
    'vs_mu', vs_mu,...
    'vs_sigma', vs_sigma,...
    't0_mu', t0_mu,...
    't0_sigma', t0_sigma);

%% Hyperparameters
% Each parameter has two hyperparms
% Build thetaprior from hyperparms structure
% pairs = reshape(fieldnames(hyperparms)', 2, numel(fieldnames(hyperparms))/2)';
% thetaprior = [];
% for i = 1:size(pairs, 1)
%     eval(sprintf('thetaprior = [thetaprior; %s, %s];', pairs{i,1}, pairs{i,2}));
% end

% Build phiprior from hyperparms structure
names = fieldnames(hyperparms);
lookForPattern = @(str,pat)~cellfun('isempty',regexp(str,pat,'once'));
abIndex = lookForPattern(names, '_sigma'); % Flag hyperparms which don't have a normal hyperprior (i.e., that have a Gamma hyperprior and therefore need an _a, _b suffix)
phiprior = [];
for i = 1:size(names, 1)
    if abIndex(i)
        eval(sprintf('phiprior = [phiprior; %s_a, %s_b];', names{i,1}, names{i,1}));
    else
        eval(sprintf('phiprior = [phiprior; %s_mu, %s_sigma];', names{i,1}, names{i,1}));
    end
end
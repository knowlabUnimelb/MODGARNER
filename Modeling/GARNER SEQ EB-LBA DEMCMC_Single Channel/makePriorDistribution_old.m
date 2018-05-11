function pd = makePriorDistribution(varname, hyperprior)

switch varname
    case {'w', 'wc', 'wp', 'b1', 'b2', 'b3', 'b4', 'b5'}
        pd = makedist('Beta', 'a', hyperprior(1), 'b', hyperprior(2)); % Sample starting points from a beta
    case {'c', 'cf', 'bh' ,'A', 'bMa1', 'bMa2', 'vs', 't0',...
            'w_a', 'w_b', 'wc_a', 'wc_b', 'c_mu', 'cf_mu',...
            'b1_a', 'b1_b', 'b2_a', 'b2_b', 'b3_a', 'b3_b', 'b4_a', 'b4_b', 'b5_a', 'b5_b',...
            'wp_a', 'wp_b', 'bh_mu', 'A_mu', 'bMa1_mu', 'bMa2_mu', 'vs_mu', 't0_mu'}
        pd = makedist('Normal', 'mu', hyperprior(1), 'sigma', hyperprior(2)); % Sample starting points from a truncated normal
        pd = truncate(pd, 0, Inf);
    case {'c_sigma', 'cf_sigma', 'bh_sigma', 'A_sigma',...
            'bMa1_sigma', 'bMa2_sigma', 'vs_sigma', 't0_sigma'}
        pd = makedist('Gamma', 'a', hyperprior(1), 'b', hyperprior(2));
end
function pdf = getPhiPriorDensity(varname, xval, hyperprior)
% See makePriorDistribution

switch varname
    case {'w', 'wc', 'b1', 'b2', 'b3', 'b4', 'b5', 'wp'}
        pdf = logbetapdf(xval, hyperprior(:,1), hyperprior(:,2));
    case {'c', 'cf', 'bh' ,'A', 'bMa1', 'bMa2', 'vs', 't0',...
            'w_a', 'w_b', 'wc_a', 'wc_b', 'c_mu', 'cf_mu',...
            'b1_a', 'b1_b', 'b2_a', 'b2_b', 'b3_a', 'b3_b', 'b4_a', 'b4_b', 'b5_a', 'b5_b',...
            'wp_a', 'wp_b', 'bh_mu', 'A_mu', 'bMa1_mu', 'bMa2_mu', 'vs_mu', 't0_mu'}
        pdf = logTruncatedNormalDistribution(xval, hyperprior(:,1), hyperprior(:,2), 0, Inf);
    case {'c_sigma', 'cf_sigma', 'bh_sigma', 'A_sigma',...
            'bMa1_sigma', 'bMa2_sigma', 'vs_sigma', 't0_sigma'}
        pdf = loggampdf(xval, hyperprior(:,1), hyperprior(:,2));
end
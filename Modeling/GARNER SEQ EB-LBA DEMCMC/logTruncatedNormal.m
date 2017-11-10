function y = logTruncatedNormal(x, mu, sigma, a, b)
    x1 = (x - mu)./sigma;
    f1 = -log(sigma) + logStandardNormalPdf(x1);
    
    if a == -inf
        f2 = 0;
    else
        a1 = (a - mu)./sigma;
        f2 = standardNormalCdf(a1);
    end
    
    if b == inf
        f3 = 1;
    else
        b1 = (b - mu)./sigma;
        f3 = standardNormalCdf(b1);
    end
    
    y = f1 - log(f2 + f3);

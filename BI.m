function [tstruct, varargout] = BI(y, x, ainu)
% Solves least square problem using Chave-Thompson bounded influence
% algorithm.
%
% Input Variables:
% y - n by 1 response vector (required)
% x - n by p predictor matrix (required)
% ainu - leverage point weighting scale factor
%
% Output Variables:
% tstruct - structure contaning number of iterations and regression
%                 parameters, se, t-statistic and p-value for each
% rstruct - structure contaning regression residuals, weights and
%                 normalized sum of squared residuals for each iteration
% pstruct - structure containing r^2, and anova, runs test and
%                 Durbin-Watson test p-values for each iteration
%
% Empircally, the leverage point weight scale factor ainu should usually be
% taken as the 0.995 quantile of the appropriate beta distribution. This is
% default. However, if the data are sustematically long tailed, a larger
% quantile may be needed.

n = length(y);
p = size(x, 2);
w = ones(size(y));

if ~exist("ainu","var")
    ainu = n/p*betainv(0.995, p, n-p); % default
end

[b, bse, t, tprob] = Treg(y, x, w);

r = y - x*b;
d = iqr(r)/1.34903;
ovar = r'*r/n;
rss = ovar;

what = ones(size(y));
wsum = n;
swhat = exp(-0.3068528*exp(-ainu^2));
hat = diag(x*inv(x'*x)*x');

bs = b;
bses = bse;
ts = t;
tprobs = tprob;

rs = [];
ws = [];
whats = [];
r2 = [];
apval = [];
rpval = [];
dwpval = [];

for k = 1:15
    wh = 1./max(1, abs(r)/(1.5*d));
    ws(:,k) = wh;
    t1 = ainu*(wsum*hat/2 - ainu);
    what = what.*exp(-0.3068528*exp(t1))/swhat;
    whats(:,k) = what;
    w = wh.*what;
    a = Anovar(y, x, w);
    r2 = [r2 a(1,2)/a(3,2)];
    apval = [apval a(1,5)];

    [b, bse, t, tprob] = Treg(y, x, w);

    bs = [bs b];
    bses = [bses bse];
    ts = [ts t];
    tprobs = [tprobs tprob];

    [h, pval] = runstest(w.*r);
    
    rpval = [rpval pval];
    wx = repmat(w, 1, p).*x;
    
    pval = dwtest(w.*r, wx);
    dwpval = [dwpval pval];
    
    hat = diag(wx*inv(wx'*wx)*wx');
    jj = find(w ~= 0);
    r = y - x*b;
    rs(:,k) = r;
    
    d = iqr(r(jj))/1.34903;
    wsum = sum(w);
    nvar = sum((w.*r).^2)/wsum;
    rss = [rss nvar];

    if abs(nvar - ovar) <= 0.01*ovar
        break
    end

    ovar = nvar;

end

tstruct.b = bs;
tstruct.bse = bses;
tstruct.t = ts;
tstruct.pval = tprobs;

if nargout > 1
    rstruct.res = rs;
    rstruct.w = ws;
    rstruct.what = whats;
    rstruct.rss = rss;
    varargout{1} = rstruct;
    varargout{2} = [];
end

if nargout > 2
    pstruct.r2 = r2;
    pstruct.anova = a;
    pstruct.runs = rpval;
    pstruct.dw = dwpval;
    varargout{2} = pstruct;
end

end
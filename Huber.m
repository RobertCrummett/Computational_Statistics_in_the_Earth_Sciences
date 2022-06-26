function [tstruct, varargout] = Huber(y, x)
% Solves least squares problem using Huber robust algorithm input variables
% 
% Input Variables:
% y - n by 1 response vector (required)
% x - n by p predictor matrix (required)
%
% Output Variables:
% tstruct - Structure containing number of iterations and regression
%                 parameters, se, t-statistic and p-value for each
% rstruct - Structure containing regression residuals, weights and
%                 normalized sum of squared residuals for each iteration
% pstruct - Structure containing r^2, and anova, runs test and
%                 Durbin-Watson test p-values for each iteration

n = length(y);
p = size(x, 2);
w = ones(size(y));

[b, bse, t, tprob] = Treg(y, x, w);

r = y - x*b;
d = iqr(r)/1.34903;
ovar = r'*r/n;
rss = ovar;
bs = b;
bses = bse;
ts = t;
tprobs = tprob;
rs = [];
ws = [];
r2 = [];
apval = [];
rpval = [];
dwpval = [];

for k = 1:15 % allow 15 iterations max
    w = 1./max(1, abs(r)/(1.5*d));
    ws(:,k) = w;
    a = Anovar(y, x, w);
    r2 = [r2 a(1,2)/a(3,2)];
    apval = [apval a(1,5)];

    [b ,bse, t, tprob] = Treg(y, x, w);
    bs = [bs b];
    bses = [bses bse];
    ts = [ts t];
    tprobs = [tprobs tprob];

    r = y - x*b;
    rs(:,k) = r;

    [~, pval] = runstest(w.*r);
    rpval = [rpval pval];
    
    pval = dwtest(w.*r, repmat(w, 1, p).*x);
    dwpval = [dwpval pval];
    d = iqr(r)/1.34903;
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
tstruct.tprob = tprobs;

if nargout > 1
    rstruct.res = rs;
    rstruct.w = ws;
    rstruct.rss = rss;
    varargout{1} = rstruct;
    varargout{2} = [];
end

if nargout > 2
    pstruct.r2 = r2;
    pstruct.anova = apval;
    pstruct.runs = rpval;
    pstruct.dw = dwpval;
    varargout{2} = pstruct;
end

end
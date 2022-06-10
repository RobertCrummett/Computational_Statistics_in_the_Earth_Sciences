function [M, sp2] = Bartlett(varargin)
    % credit to Wikipedia
    % this get the chisqstat exactly, unlike Dr. Chave's function
    % R N Crummett 6/6/2022
    k = nargin;
    n = zeros(1,k);
    s = zeros(1,k);

    for i = 1:k
        n(i) = length(varargin{i});
        s(i) = std(varargin{i});
    end
    
    N = sum(n);

    sp2 = (1/(N - k))*sum((n - 1).*(s.^2));

    M = ((N-k)*log(sp2) - sum((n -1).*log(s.^2)))./(1 + (1/(3*(k - 1)))*(sum(1./(n - 1)) - 1/(N - k)));
end

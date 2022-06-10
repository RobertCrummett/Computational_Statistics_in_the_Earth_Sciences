function [Sum] = Hypergeometric2f1(a, b, c, x)
%Computes Gauss hypergeometric function using Pade approximations to
%accelerate convergence

rerr = 1e-8; %relative error
aerr = 1e-35; %absolute error
last = 0;
s = [];

if (a < 0) && ((a - fix(a)) == 0)
    kup = -a;
elseif (b < 0) && ((b - fix(b)) == 0)
    kup = -b;
else
    kup = 1000; %upper limit is arbitrary and should not be reached
end

if x ~= 1
    for k = 0:kup
        pocha = Pochhammer(a, k);
        pochb = Pochhammer(b, k);
        pochc = Pochhammer(c, k);
        s = [s pocha*pochb*x^k/(pochc*gamma(k+1))];
        Sum = Padesum(s);
        if abs (Sum - last) <= rerr*abs(Sum) + aerr
            return
        end
        last = Sum;
    end
else
    Sum = gamma(c)*gamma(c-a-b)/(gamma(c-a)*gamma(c-b));
end
end

function [Result] = Pochhammer(x, k)
%computes Pochhammer's symbol
if k == 0
    Result = 1;
    return
else
    i = 0:k -1;
    s = x + i;
    Result = prod(s);
    return
end
end

function [Cf] = Padesum(s)
%Computes sum from 1 to n of s(i) using Pade approximant implemented with
%continued fraction expansion; see Z. Naturforschung, 33a, 402-417, 1978.
%
%Input variable
%s--series of values to be summed. may be complex
%Output variable
%Cf--sum of the series

n = length(s);
D = [];
x = zeros(1,n);
d = [];
t = zeros(1,n);
D(1) = s(1);
d(1) = D(1);
if n == 1
    Cf = d(1);
    return
end

D(2) = s(2);
d(2) = -D(2)/D(1);
if n == 2
    Cf = d(1)/(1 + d(2));
    return
end

for i = 3:n
    L = 2*fix((i-1)/2);
    %update x vector
    x(L:-2:4) = x(L-1:-2:3) + d(i - 1)*x(L-2:-2:2);
    x(2) = x(1) + d(i-1);
    %interchange odd and even parts
    t(1:2:L-1) = x(1:2:L-1);
    x(1:2:L-1) = x(2:2:L);
    x(2:2:L) = t(1:2:L-1);
    %compute cf coefficent
    D(i) = s(i) + s(i - 1:-1:i - L/2)*(x(1:2:L-1)).';
    d(i) = -D(i)/D(i-1);
end
%evaluate continued fraction
Cf = 1;
for k = n:-1:2
    Cf = 1 + d(k)/Cf;
end
Cf = d(1)/Cf;
end

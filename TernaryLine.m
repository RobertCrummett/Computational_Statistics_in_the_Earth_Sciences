function TernaryLine(comp, varargin)
%Plots a compositional line given the 3 part composition in comp
xo = (comp(1,1) + 2*comp(1,2))/(2*sum(comp(1,:)));
yo = sqrt(3)*comp(1,1)/(2*sum(comp(1,:)));

if nargin == 1
    col = 'k';
    lw = 1;
elseif nargin == 2
    col = varargin{1};
    lw = 1;
elseif nargin == 3
    col = varargin{1};
    lw = varargin{2};
end

for i=2:length(comp)
    x = (comp(i,1) + 2*comp(i,2))/(2*sum(comp(i,:)));
    y = sqrt(3)*comp(i,1)/(2*sum(comp(i,:)));
    plot([xo x],[yo y], col, 'linewidth', lw)
    xo=x;
    yo=y;
end
end


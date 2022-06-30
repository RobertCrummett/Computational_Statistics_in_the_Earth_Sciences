function TernaryPlot(comp, varargin)
%Plots the compositional data in comp using blue + by default
%varargin{1} specifies alternate symbols (e.g.,'d') and varargin{2}
%specifies alternate collors (e.g.,'k')
xx = Comp_to_Cart(comp);
if nargin == 1
    mrk = '+';
    mcol = 'b';
elseif nargin == 2
    mrk = varargin{1};
    mcol = 'b';
else
    mrk = varargin{1};
    mcol = varargin{2};
end
for i = 1:length(xx)
    plot(xx(i,1), xx(i,2), 'Marker', mrk, 'MarkerSize', 8, ...
        'MarkerEdgeColor', mcol, 'MarkerFaceColor', mcol)
end
end
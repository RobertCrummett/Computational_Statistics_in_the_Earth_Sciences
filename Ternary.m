function [h] = Ternary(varargin)
%sets up axes for ternary diagram
%arguments
%ntick - number of tick marks on the axes (default is 11); if set to zero,
%then tick values are contained in varargin{3}
%composition - composition by which the gridlines are perturbed
if nargin == 0
    ntick = 11;
else
    ntick = varargin{1};
end
hold off
h = fill([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'w', 'linewidth', 2);
axis image
axis off
hold on
label = linspace(0, 1, ntick);
if nargin < 2
    label = linspace(0, 1, ntick);
    for i = 1:ntick
        plot([label(i) 0.5*(1 + label(i))], [0 0.866*(1 - label(i))], ...
             ':k', 'linewidth', 0.75)
        plot([0.5*(1 - label(i)) 1 - label(i)], ...
            [0.866*(1 - label(i)) 0], ':k', 'linewidth', 0.75)
        plot([0.5*label(i) 1 - 0.5*label(i)], ...
             [0.866*label(i) 0.866*label(i)], ':k', 'linewidth', 0.75)
        text(1 - 0.5*label(i) + .025, sqrt(3)*label(i)/2 + 0.01, ...
             num2str(label(i), 3), 'FontSize', 14, 'FontWeight', 'bold');
        str = num2str(label(i), 3);
        text(label(i) - 0.02*length(str)/2, -.025, str, 'FontSize', 14,...
             'FontWeight', 'bold');
        str = num2str(label(i), 3);
        text(0.5*(1 - label(i)) - 0.015*(length(str) - 1) -.035, ...
             sqrt(3)/2*(1 - label(i)) + .01, str, 'FontSize', 14, ...
             'FontWeight', 'bold');
    end
elseif nargin == 2
    label = linspace(0, 1, ntick);
    y = varargin{2};
    for i = 1:ntick
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([label(i) 0]), y));
         
         x2 = Comp_to_Cart(Perturbation( ... 
              Cart_to_Comp([0.5*(1 +label(i)) 0.866*(1 - label(i))]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([0.5*(1 - label(i)) 0.866*(1 - label(i))]), y));
         x2 = Comp_to_Cart(Perturbation( ...
             Cart_to_Comp([1 - label(i) 0]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([0.5*label(i) ...
              0.866*label(i)]), y));
         x2 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([1 - 0.5*label(i) 0.866*label(i)]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([1 - 0.5*label(i) 0.866*label(i)]), y));
         text(x1(1) + .025, x1(2) + 0.01, num2str(label(i), 3), ...
              'FontSize', 14, 'FontWeight', 'bold');
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([label(i) 0]), y));
         str = num2str(label(i), 3);
         text(x1(1) - 0.02*length(str)/2, x1(2)-.025, str, 'FontSize', ...
             14, 'FontWeight', 'bold');
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([0.5*(1 - label(i)) 0.866*(1 - label(i))]), y));
         str = num2str(label(i), 3);
         text(x1(1) - .015*(length(str) - 1) - .035, x1(2) + .01, str, ...
             'FontSize', 14, 'FontWeight', 'bold');
    end
else
    label = varargin{3};
    ntick = length(label);
    for i = 1:ntick
         y = varargin{2};
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([label(i) 0]), y));
         x2 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([0.5*(1 +label(i)) 0.866*(1 - label(i))]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([0.5*(1 - label(i)) 0.866*(1 - label(i))]), y));
         x2 = Comp_to_Cart(Perturbation( ...
             Cart_to_Comp([1 - label(i) 0]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([0.5*label(i) ...
              0.866*label(i)]), y));
         x2 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([1 - 0.5*label(i) 0.866*label(i)]), y));
         plot([x1(1) x2(1)], [x1(2) x2(2)], ':k', 'linewidth', 0.75)
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([1 - 0.5*label(i) 0.866*label(i)]), y));
         text(x1(1) + .025, x1(2) + 0.01, num2str(label(i), 3), ...
              'FontSize', 14, 'FontWeight', 'bold');
         x1 = Comp_to_Cart(Perturbation(Cart_to_Comp([label(i) 0]), y));
         str = num2str(label(i), 3);
         text(x1(1) - 0.02*length(str)/2, x1(2)-.025, str, 'FontSize', ...
             14, 'FontWeight', 'bold');
         x1 = Comp_to_Cart(Perturbation( ...
              Cart_to_Comp([0.5*(1 - label(i)) 0.866*(1 - label(i))]), y));
         str = num2str(label(i), 3);
         text(x1(1) - .015*(length(str) - 1) - .035, x1(2) + .01, str, ...
             'FontSize', 14, 'FontWeight', 'bold');
    end
end
end

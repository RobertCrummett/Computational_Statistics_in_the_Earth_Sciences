function [Result] = Ilr(x, Psi)
% Computes the ilr of the composition x using the contrast matrix Psi
Result = Clr(x)*Psi';
end
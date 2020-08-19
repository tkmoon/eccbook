function [x,pg] = plotgauss(m,s2)
% function [x,pg] = plotgauss(m,s2)
% Plot a Gaussian function with mean m and variance s2

m1 = m-5*s2;
m2 = m+5*s2;
x = linspace(m1,m2,200);

pg = 1/sqrt(2*pi*s2) * exp(-.5*(x - m).^2 / s2);

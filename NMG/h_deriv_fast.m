function [h1, h2, h3, h12, h13, h23, h33] = h_deriv_fast(a, h10, h20, h30,...
    h120, h130, h230, h330)
% compute h1, h2, h3, h12, h13, h23, h33

a_sq = a^2;
h1 = a_sq*h10;
h2 = a_sq*h20;
h3 = a_sq*h30;
h12 = a_sq*h120;
h13 = a_sq*h130;
h23 = a_sq*h230;
h33 = a_sq*h330;

end
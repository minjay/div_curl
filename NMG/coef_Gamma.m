function [Gamma1, Gamma2] = coef_Gamma(ai, aj, bi, bj, h1, h2, h3, h12, h13, h23, h33)

Gamma1 = 1/4*(ai*aj*h1*h2+bi*aj*h2*h3-ai*bj*h1*h3-bi*bj*h3^2);
Gamma2 = -1/2*(ai*aj*h12+bi*aj*h23-ai*bj*h13-bi*bj*h33);

end
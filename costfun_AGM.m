function [E, G] = costfun_AGM(param, X, Y, K, U, lambda, sigma2, is_grad,r)
% ##=======================================================================
%  Point Set Registration with Asymmetric Gaussian Mixtures (AGMReg) Demo 
% ##=======================================================================
% 
% Author:       Gang Wang
% Date:         03/15/2015
% Email:    gwang.cv@gmail.com
%
%--------------------------------------------------------------------------
[N, D] = size(X);
M = size(U, 2);
C = reshape(param, [M D]);
E0 = 1/( 2^D * pi^(D/2) * (sigma2*(((r+1)/2)^2))^(D/4) );
E1 = E0 + lambda/2 * trace(C'*K*C);
a = -2 / N / (2*pi*sigma2)^(D/2) / ((r+1)/2)^D;
V = U*C - Y;  
Va=sum(V,2);
F0 = exp(-sum(V.^2, 2) / (2*sigma2));
F1 = exp(-sum(V.^2, 2) / (2*sigma2*r*r));
F=((Va<=0).*F0+(Va>0).*F1); 
E = E1 + a * sum(F);
G = [];
if is_grad
    Vs=repmat(sum(V,2),[1 D]);    
    FF=repmat(F, [1 D]);
    Vb= V .* FF;
    Vc=(Vs<=0).*Vb/sigma2+(Vs>0).*Vb/(sigma2*r*r);
    G = - a * U' * ( Vc ) + lambda * K * C;
    G = G(:);
end

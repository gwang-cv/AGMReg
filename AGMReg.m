function [Y_reg]=AGMReg(X,Y,iterNum)
% ##=======================================================================
%  Point Set Registration with Asymmetric Gaussian Mixtures (AGMReg) Demo 
% ##=======================================================================
% 
% Author:       Gang Wang
% Date:         03/15/2015
% Email:    gwang.cv@gmail.com
%
%--------------------------------------------------------------------------
% initialize
normalize = 1;
visualize =1;
s=1;
[n1,D]=size(X);
[n2,D]=size(Y);
if visualize
    h1=figure(1);
    set(gca,'FontSize',12,'FontName','Times','FontWeight','bold');
    if D==2
        plot(X(1:n1,1),X(1:n1,2),'b+',Y(:,1),Y(:,2),'ro','MarkerSize', 5, 'LineWidth',1);
    else
        plot3(X(1:n1,1),X(1:n1,2),X(1:n1,3),'b+',Y(:,1),Y(:,2),Y(:,3),'ro','MarkerSize', 5, 'LineWidth',1);
    end
    title('Initialization');
end
T =0.5;
it_total=1;
anneal_rate=0.93;
if nargin<3
    T_finalfac=500;
else
    T_finalfac = iterNum;
end
T_final=T/T_finalfac; 
R=1.5;
lamda1_init = 1;
lamda2_init = 0.01;
[xmax, dim] = size(X); X = X (1:1:xmax, :); [xmax, dim] = size(X);
[ymax, tmp] = size(Y); Y = Y (1:1:ymax, :); [ymax, tmp] = size(Y);
z = X;
c_tps = zeros (xmax,dim+1);
d_tps = eye   (dim+1, dim+1);
w     = zeros (xmax+dim+1, dim);
m              = ones (xmax, ymax) ./ (xmax * ymax);
T0             = max(X(:,1))^2;
moutlier       = 1/sqrt(T0)*exp(-1);      
m_outliers_row = ones (1,ymax) * moutlier;
m_outliers_col = ones (xmax,1) * moutlier;
Xk=X;
Yk=Y;
figure(2);
while s
    perT_maxit=5;
    for i=1:perT_maxit  
        m = CalCorres (Xk, Y, T, 'mix-rpm', m_outliers_row, m_outliers_col,it_total);
        vy = m * Y ./ ( (sum(m'))' * ones(1,dim));
        lamda1 = lamda1_init*length(X)*T;
        lamda2 = lamda2_init*length(X)*T;
        [c_tps, d_tps, w] = cMIX_calc_transformation ('tps',lamda1, lamda2, 1, X, vy, z);
        [vx] = cMIX_warp_pts ('tps', X, z, c_tps, d_tps, w, 1);
        Xk=vx;
    end
    V=vx;
        if visualize           
            if D==2
                plot(V(1:n1,1),V(1:n1,2),'b+',Y(:,1),Y(:,2),'ro','MarkerSize', 5, 'LineWidth',1);
            else
                plot3(V(1:n1,1),V(1:n1,2),V(1:n1,3),'b+',Y(:,1),Y(:,2),Y(:,3),'ro','MarkerSize', 5, 'LineWidth',1);
            end
             drawnow           
        end
    T = T * anneal_rate;
    if T < T_final; s = 0; end;
    it_total = it_total + 1;
end  
for ind=1:10
    perT_maxit=5;
    for i=1:perT_maxit   
        m = CalCorres (Xk, Y, T, 'mix-rpm', m_outliers_row, m_outliers_col,it_total);
        vy = m * Y ./ ( (sum(m'))' * ones(1,dim)); 
        lamda1 = lamda1_init*length(X)*T;
        lamda2 = lamda2_init*length(X)*T;
        [c_tps, d_tps, w] = cMIX_calc_transformation ('tps',lamda1, lamda2, 1, X, vy, z);
        [vx] = cMIX_warp_pts ('tps', X, z, c_tps, d_tps, w, 1);
        Xk=vx;
    end
    X2 = Xk;
    Y2=vy;
    normal.xm=0; normal.ym=0;
    normal.xscale=1; normal.yscale=1;
    if normalize
        [nX, nY, normal]=norm2s(X2,Y2);
        [idt, V, T,R] = AGMcost(nX, nY-nX, 0.5, R,T); 
    end
    if normalize
        V=(V+nX)*normal.yscale+repmat(normal.ym,size(Y2,1),1);
    end
    if visualize
        if D==2
            plot(V(1:n1,1),V(1:n1,2),'b+',Y(:,1),Y(:,2),'ro','MarkerSize', 5, 'LineWidth',1);
        else
            plot3(V(1:n1,1),V(1:n1,2),V(1:n1,3),'b+',Y(:,1),Y(:,2),Y(:,3),'ro','MarkerSize', 5, 'LineWidth',1);
        end
        drawnow
    end
    Xk = V;
end
title('After registration');
Y_reg=V;
end

%% Funtions
function [idt, V , sigma2,r] = AGMcost(X, Y, thresh,r,sigma2)
[N1, D] = size(X);
[N2, D] = size(Y);
N=min(N1,N2);
is_grad = 1;
beta =0.8; 
lambda = 0.1; 
arfa=0.93;
iter=2;
if N>=15
    n_ker = 15;
else
    n_ker=N;
end
[Q,K]=lowrankQS(X, beta, n_ker, 1);
U=Q*K;
x0 = zeros(n_ker*D, 1);
options = optimset( 'display','off', 'MaxIter', 100);
if is_grad
    options = optimset(options, 'GradObj', 'on');
end
param = fminunc(@(x)costfun_AGM(x, X, Y, K, U, lambda, sigma2, is_grad,r), x0, options); 
for ii = 1:iter
    sigma2 = sigma2*arfa;
    r=r*0.85;
    param = fminunc(@(x)costfun_AGM(x, X, Y, K, U, lambda, sigma2, is_grad,r), param, options);
end
C = reshape(param, [n_ker D]);
V=U*C;
if N1~=N2
    V0 = Y- [V;zeros(abs(N2-N1),2)] ;
else
    V0 = Y- V;
end
Va=sum(V0,2);
F0 = exp(-sum(V0.^2, 2) / (2*sigma2));
F1 = exp(-sum(V0.^2, 2) / (2*sigma2*r*r));
Pb=((Va>0).*F0+(Va<=0).*F1);
idt = find(Pb > thresh);
end

function [m, m_outliers_row, m_outliers_col] = CalCorres(vx, y, T, m_method, m_outliers_row, m_outliers_col, it_total, icp_sigma)
[xmax,dim] = size(vx);
[ymax,dim] = size(y);
if strcmp (m_method, 'mix-rpm')    
    y_tmp = zeros (xmax, ymax);
    for it_dim=1:dim
        y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
    end;
    m_tmp = 1/sqrt(T) .* exp (-y_tmp/T);
    m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001;
    m = m_tmp;
    [m, junk1, junk2] = cMIX_normalize_m (m_tmp, m_outliers_col, m_outliers_row);
end
end

function [vx] = cMIX_warp_pts (trans_type, x, z, c_tps, d_tps, w, sigma_kernel);

switch (trans_type)
  case 'tps'
    vx = ctps_warp_pts (x, z, c_tps, d_tps); 
  case 'rbf'
    vx = crbf_warp_pts (x, z, w, sigma_kernel);
  otherwise; 
    disp ('# ERROR #: cMIX_warp_pts -- wrong input!');
end
end

function [m, m_outlier_col, m_outlier_row] = cMIX_normalize_m (m, ...
    m_outlier_col, m_outlier_row);

norm_threshold = 0.05;
norm_maxit     = 10;

[xmax, ymax] = size(m);
  
norm_it = 0;
while (1 > 0)

  sx = sum(m')' + m_outlier_col;
  m  = m ./ (sx * ones(1,ymax));
  m_outlier_col = m_outlier_col ./sx;
  
  sy = sum(m) + m_outlier_row;
  m  = m ./ (ones(xmax,1)*sy);
  m_outlier_row = m_outlier_row ./sy;
  err = ((sx-1)'*(sx-1) + (sy-1)*(sy-1)')/(length(sx)+length(sy));
  if err < (norm_threshold .* norm_threshold); break, end
  norm_it = norm_it + 1;
  if norm_it >= norm_maxit; break; end
  
end
end

function [c_tps, d_tps, w] = cMIX_calc_transformation (transformation_type, ...
    lamda1, lamda2, sigma_kernel, x, vy, z);

c_tps = [];
d_tps = [];
w     = [];

switch (transformation_type)
  case 'tps'
    [c_tps,d_tps]  = ctps_gen (x, vy, lamda1, lamda2);
  case 'rbf'
    [phi, w] = crbf_gen (x, vy, z, lamda1, lamda2, sigma_kernel);
  otherwise; 
    disp ('# ERROR #: cMIX_calc_transformation -- wrong input!');
end
end
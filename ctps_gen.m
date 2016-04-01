% -------------------------------------------------------------------
% ctps_gen.m
% -------------------------------------------------------------------
% Purpose: Generate all parameters for TPS.
%
% Usage: 
% [K]             = ctps_gen (x);
% [K]             = ctps_gen (x, y);
% [c,d]           = ctps_gen (x, y, lamda1);
% [c,d]           = ctps_gen (x, y, lamda1, lamda2);
% [q1,q2,R,K]     = ctps_gen (x);
% [q1,q2,R,K,c,d] = ctps_gen (x, y, lamda1);

function [o1,o2,o3,o4,o5,o6] = ctps_gen (x,y,lamda1,lamda2);

% check input.
if nargin <= 1 | nargin >= 5
  disp ('# ERROR #: ctps_gen -- wrong input !');
  help ctps_gen; return;
end;

o1 = [];
o2 = [];
o3 = [];
o4 = [];
o5 = [];
o6 = [];

% --- [K] = ctps_gen (x) --------------------------------------------
if nargin == 1 & nargout == 1
  [n, dim] = size (x); x = [ones(n,1), x];
  [K] = ctps_gen_K (x,x); 
  
  o1 = K; 

% --- [K] = ctps_gen (x,y) ------------------------------------------
elseif nargin == 2 & nargout == 1
  [n, dim] = size (x); x = [ones(n,1), x];
  [m, dim] = size (y); y = [ones(m,1), y];
  [K] = ctps_gen_K (x,y); 
  
  o1 = K; 

% --- [c,d] = ctps_gen (lamda, t, y) --------------------------------
elseif nargin == 3 & nargout == 2
  [n, dim] = size (x); x = [ones(n,1), x];
  [m, dim] = size (y); y = [ones(n,1), y];

  [K]       = ctps_gen_K (x,x);
  [q1,q2,R] = ctps_gen_qr     (x);
  
  [c,d]     = ctps_gen_cd     (lamda1,q1,q2,R,K,y);

  o1 = c;
  o2 = d;

% --- [c,d] = ctps_gen (x, y, lamda1, lamda2) ----------------------
elseif nargin == 4 & nargout == 2
  [n, dim] = size (x); x = [ones(n,1), x];
  [m, dim] = size (y); y = [ones(n,1), y];

  [K]       = ctps_gen_K (x, x);
  [q1,q2,R] = ctps_gen_qr     (x);
  
  [c,d]     = ctps_gen_cd_regularized (lamda1,lamda2,q1,q2,R,K,y);

  o1 = c;
  o2 = d;

% --- [q1,q2,R,K] = ctps_gen (x) -----------------------------------
elseif nargin == 1 & nargout == 4
  [n, dim] = size (x); x = [ones(n,1), x];

  [K]       = ctps_gen_K (x, x);
  [q1,q2,R] = ctps_gen_qr     (x);
  
  o1 = q1;
  o2 = q2;
  o3 = R;
  o4 = K;

% --- [q1,q2,R,K,c,d] = ctps_gen (x, y, lamda1) --------------------
elseif nargin == 3 & nargout == 6
  [n, dim] = size (x); x = [ones(n,1), x];
  [m, dim] = size (y); y = [ones(n,1), y];

  [K]       = ctps_gen_K (x, x);
  [q1,q2,R] = ctps_gen_qr     (x);
  
  [c,d]     = ctps_gen_cd (lamda1,q1,q2,R,K,y);

  o1 = q1;
  o2 = q2;
  o3 = R;
  o4 = K;
  o5 = c;
  o6 = d;
else
  disp ('# ERROR #: ctps_gen -- wrong input!');
  help ctps_gen;
end;

function [q1,q2,R] = ctps_gen_qr (x);

[n,M] = size (x);

[q,r]   = qr(x);
q1      = q(:, 1:M);
q2      = q(:, M+1:n);
R       = r(1:M,1:M);


function [K] = ctps_gen_K (x,z);

% Format:
[n, M] = size (x); 
[m, M] = size (z);
dim    = M  - 1;

% calc. the K matrix.
% 2D: K = r^2 * log r
% 3D: K = -r
K= zeros (n,m);

for it_dim=1:dim
  tmp = x(:,it_dim+1) * ones(1,m) - ones(n,1) * z(:,it_dim+1)';
  tmp = tmp .* tmp;
  K = K + tmp;
end;
  
if dim == 2
  mask = K < 1e-10; % to avoid singularity.
  K = 0.5 * K .* log(K + mask) .* (K>1e-10);
else
  K = - sqrt(K);
end;
  



function [c,d] = ctps_gen_cd (lamda1,q1,q2,R,K,y);

[n,M] = size(y);

gamma = inv (q2'*K*q2 + lamda1*eye(n-M, n-M)) * q2' * y;
c     = q2 * gamma;
d     = inv(R) * q1' * (y-K*q2*gamma);



function [c, d] = ctps_gen_cd_regularized (lamda1,lamda2,q1,q2,R,K,y);

[n,M] = size(y);
dim = M - 1;

gamma =  inv (q2'*K*q2 + lamda1*eye(n-M)) * q2' * y;
c     = q2*gamma;

% add regularization for "d" as well:
% d = inv(R) * q1' * (y-K*q2*gamma);
A = inv(R'*R + lamda2 * eye(length(R),length(R))) * ( R'*q1'*(y-K*q2*gamma) - R'*R);
d = A + eye(dim+1,dim+1);






function  [X, Y, normal] =norm2s(x,y)
%
% Date:         03/15/2015
% Email:    gwang.cv@gmail.com
%
% norm2s nomalizes the data so that it has zero means and unit covariance
%--------------------------------------------------------------------------

x = double(x);
y = double(y);

n=size(x,1);
m=size(y,1);

normal.xm=mean(x);
normal.ym=mean(y);

x=x-repmat(normal.xm,n,1);
y=y-repmat(normal.ym,m,1);

normal.xscale=sqrt(sum(sum(x.^2,2))/n);
normal.yscale=sqrt(sum(sum(y.^2,2))/m);

X=x/normal.xscale;
Y=y/normal.yscale;
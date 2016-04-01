function  DisplayPoints_wg(Model, Scene, dim)
% 
% Date:         03/15/2015
% Email:    gwang.cv@gmail.com
%
%--------------------------------------------------------------------------

if dim==2
set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');
plot(Scene(:,1),Scene(:,2),'ro','MarkerSize', 10, 'LineWidth',2);
hold on;
plot(Model(:,1),Model(:,2),'b+','MarkerSize', 10, 'LineWidth',2);
axis equal;
end

if dim==3
set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');
X=Scene;
Y=Model;
plot3(X(:,1),X(:,2),X(:,3),'ro',Y(:,1),Y(:,2),Y(:,3),'b+','MarkerSize', 10, 'LineWidth',2);
axis equal;
end
pbaspect([1,1,1]);
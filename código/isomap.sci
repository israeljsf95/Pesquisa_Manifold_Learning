n = 1000; 
x = rand(2,n);
v = 3*%pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(2,:) = 20*x(2,:);
X(1,:) = -cos(v).*v;
X(3,:) = sin(v).*v;
ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
clf;
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled');

view(v1,v2);axis('equal');axis('off');

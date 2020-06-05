clear, close all
n = 350;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;


ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
figure
subplot(311)
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

% matriz de distancias euclidianas
D = euclid_dist_linalg(X',2);
%Escolhendo os k vizinhos
K = 7;
%Organiza D 
G = zeros(size(D));
for i =1:n
   xaux = D(i,:);
   [~, I] = sort(xaux);
   for k = 2:K+1
       G(i,I(k)) = xaux(I(k));
   end
end
D_g = zeros(size(G));
D_g(G==0) = Inf;
D_g = G + D_g;
%Algoritmo de Floyd
for i=1:n
    Aux = repmat(D_g(:,i), 1, n)+repmat(D_g(i,:), n, 1);
    D_g = min(D_g,Aux); 
end

D_g(D_g==Inf) = 0;
D_g = (D_g + D_g')/2;
M = eye(size(D_g)).*(diag(diag(D_g)));
D_g2 = D_g - M;
A = (-1/2)*D_g2.^2;  
H2 = eye(n) - (1/n)*(ones(n,n));
B2 = (H2*A*H2);



eta = 0.001;
d = 2;
p = rand(d, size(B2,1));
cont = 2;
while (1)
    MM = euclid_dist_linalg(p', 2);%matriz de distancias dos pontos
    for k =1:size(B2,1)
%         k = randi([1,n]);
        nab_J = zeros(d,1);
        v = MM(k,:)';
        %Pegando as distancias que não sofrem o processo de Centralização
        delta = D_g2(k,:)';
        w = v - delta;
        norm_w = norm(w)*(1.3);
        for i = 1:d
            for j = 1:size(B2,2)
               if (j==k)
                   continue
               end
               nab_J(i) = nab_J(i)+ (p(i,k)-p(i,j))*(1-(delta(j)/v(j))); 
            end 
        end
        nab_J = nab_J*1/(norm_w);
        
        p(:,k) = p(:,k) - eta*nab_J;

%         subplot(312)
%         drawnow
%         plot3(p(1,:),p(2,:), p(3,:),'.')
%         axis square
    end
%     erro(cont) = sum((D(:)-M(:)).^2)/size(D,1);
%Esta parte mostra no console o valor de J diminuindo
   J(cont) = norm_w;
   disp(J(end))
   cont = cont + 1;
%    subplot(313)
%    drawnow
%    plot(J)
end

% scatter3(p(1,:), p(2,:), p(3,:), ms, v, 'filled')
% colormap jet(256);
scatter(p(1,:), p(2,:), ms, v, 'filled')
colormap jet(256);
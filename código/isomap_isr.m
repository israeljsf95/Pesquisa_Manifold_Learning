clear, close all
n = 800;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;


% 
% X2(1,:) = sin(2*pi*440*linspace(0,0.125,n));
% X2(2,:) = sin(2*pi*554.3653*linspace(0,0.125,n));
% X2(3,:) = sin(2*pi*659.2551*linspace(0,0.125,n));
% 
% X3(1,:) = sin(2*pi*do*linspace(0,0.125,n));
% X3(2,:) = sin(2*pi*mi*linspace(0,0.125,n));
% X3(3,:) = sin(2*pi*sol*linspace(0,0.125,n));
% 
% X4(1,:) = sin(2*pi*do*linspace(0,0.125,n));
% X4(2,:) = sin(2*pi*mib*linspace(0,0.125,n));
% X4(3,:) = sin(2*pi*sol*linspace(0,0.125,n));
% 
% X5(1,:) = sin(2*pi*440*linspace(0,0.125,n));
% X5(2,:) = sin(2*pi*do*linspace(0,0.125,n));
% X5(3,:) = sin(2*pi*mi*linspace(0,0.125,n));
% 

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
%aplicando o mds classico no isomap
M = eye(size(D_g)).*(diag(diag(D_g)));
D_g2 = D_g - M;
A = (-1/2)*D_g2.^2;  
H2 = eye(n) - (1/n)*(ones(n,n));
B2 = (H2*A*H2);
n_auto_val = 3;
[auto_vet, auto_val] = eigs(B2, n_auto_val, 'LR');
% [aut_val_desc, ind_aut_val] = sort(auto_val, 'descend');
s = struct();
s2 = struct();
euc_Y = struct();
for i=1:n_auto_val
    a_vet = auto_vet(1:end,1:i);
    a_val = auto_val(1:i,1:i);
%     a_vet = auto_vet(1:end,end:end-i+1);
%     a_val = auto_val(end:end-i+1,end:end-i+1);
    s(i).Y = a_vet*sqrt(a_val);
    euc_Y(i).Y = euclid_dist_linalg(s(i).Y,2);
    %Correlação entre a matriz de Distancias geodesicas, e 
    %a matriz de Distancias dos pontos reconstruidos
    J2(i) = 1 - corr2(euc_Y(i).Y, D_g2)^2;%Variância Residual
end

subplot(312)
scatter(s(2).Y(:,1), s(2).Y(:,2),ms, v, 'filled')
colormap jet(256);
subplot(313)

scatter3(s(3).Y(:,1), s(3).Y(:,2), s(3).Y(:,3), ms, v, 'filled')
colormap jet(256);

figure
grid on
plot([1:n_auto_val], J2, '-o')
title('ISOMAP Reconstrução - Swiss Roll')
xlabel('auto-valor (Dimensão)')
ylabel('Variancia Residual')

%Israel Jesus Santos Filho
%LLE
clear, close all

%%Construindo o Conjunto de Pontos
n = 1000;  
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
subplot(211)
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

%Escolha do N�mero de Vizinhos;


d = 2;
K = 7;
% matriz de distancias euclidiana entre os pontos
D = euclid_dist_linalg(X',2);
%Organiza D -> para conseguir pegar construir a matriz de pesos
%que relacionar� a reconstru��o de cada conjunto de pontos ao
%seu vizinho.
[D_ordenado, Ind] = sort(D);


%Diferentemente do ISOMAP, n�s construirremos uma matriz auxiliar
%que conter� a informa��o dos vizinhos, por�m n�o com o objetivo 
%de manter a dist�ncia, mas com objetivo de saber quais ser�o os
%os pontos na matriz de distancias que receber�o os pesos

Eta = zeros(size(D,2),K);
for i = 1:n
   Eta(i,1:K) = Ind(2:K+1,i)';
end


%matriz contendo os vizinhos

clear Aux
%construindo a matriz de pesos
%ajustar os pesos
G = zeros(size(X,2), size(X,2));
W = zeros(size(X,2),size(X,2));

%Vamos criar a matriz de Corela��o para estima��o dos pesos de re
%constru��o do mapa, a partir dai, come�a o processo de abertura
%para estima��o do mapa de redimensionamento
%Minha matriz de correla��o estimada ser� de dimens�o KxK, em que
%K � a quantidade de vizinhos
alfa = 0;
beta = 0;
for i=1:n
   x_i = X(:,i);%pegando o valor do ponto
   viz_x_i = X(:,Eta(i,:));
   C = (viz_x_i'*viz_x_i)/(n-1);%estimando a matriz de correla��o
   C_inv = inv(C + (1e6)*eye(size(C,1)));%calculando a inversa da matriz
   %Calculando o multiplicador de lagrange
   alfa = 1 - sum(C_inv*(x_i'*viz_x_i)');
   beta = sum(sum(C_inv));
   lambda = alfa/beta;%multiplicador de lagrange
   W(i,Eta(i,:)) = (C_inv*(x_i'*viz_x_i + lambda)')';
end
    

%Garantindo que os pesos somem 1
for i=1:size(W,1)
    W(i,:) = W(i,:)/sum(W(i,:));
end

%Reconstruindo o MDS cl�ssico
M = (eye(size(W,1)) - W)'*(eye(size(W,1)) - W);
% [auto_vet, auto_val] = eigs(M, 2+1, 0);%dos menores
[Y, aut_val] = eigs(M, d+1, 0);
%autovalores, e autovetores associados, aos maiores
% [auto_vet, auto_val] = eig(M);
Y = Y(:,2:d+1)'*sqrt(n)
% s = struct();
% s2 = struct();
% euc_Y = struct();
% for i=K-1:-1:1
%     a_vet = auto_vet(1:end,1:i);
% %     a_val = auto_val(1:i,1:i);
% %     a_vet = auto_vet(1:end,end:end-i+1);
% %     a_val = auto_val(end:end-i+1,end:end-i+1);
%     s(i).Y = a_vet;
% %     euc_Y(i).Y = euclid_dist_linalg(s(i).Y,2);
% end



subplot(212)
scatter(Y(1,:), Y(2,:), ms, v, 'filled')
% scatter3(Y(1,:), Y(2,:),Y(3,:), ms, v, 'filled')
colormap jet(256);


































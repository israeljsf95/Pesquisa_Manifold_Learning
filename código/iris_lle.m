%brincando com o dado iris
%visualizando com o isomap 
%
clear all, close all

iris = [5.1,3.5,1.4,0.2,1
4.9,3.0,1.4,0.2,1
4.7,3.2,1.3,0.2,1
4.6,3.1,1.5,0.2,1
5.0,3.6,1.4,0.2,1
5.4,3.9,1.7,0.4,1
4.6,3.4,1.4,0.3,1
5.0,3.4,1.5,0.2,1
4.4,2.9,1.4,0.2,1
4.9,3.1,1.5,0.1,1
5.4,3.7,1.5,0.2,1
4.8,3.4,1.6,0.2,1
4.8,3.0,1.4,0.1,1
4.3,3.0,1.1,0.1,1
5.8,4.0,1.2,0.2,1
5.7,4.4,1.5,0.4,1
5.4,3.9,1.3,0.4,1
5.1,3.5,1.4,0.3,1
5.7,3.8,1.7,0.3,1
5.1,3.8,1.5,0.3,1
5.4,3.4,1.7,0.2,1
5.1,3.7,1.5,0.4,1
4.6,3.6,1.0,0.2,1
5.1,3.3,1.7,0.5,1
4.8,3.4,1.9,0.2,1
5.0,3.0,1.6,0.2,1
5.0,3.4,1.6,0.4,1
5.2,3.5,1.5,0.2,1
5.2,3.4,1.4,0.2,1
4.7,3.2,1.6,0.2,1
4.8,3.1,1.6,0.2,1
5.4,3.4,1.5,0.4,1
5.2,4.1,1.5,0.1,1
5.5,4.2,1.4,0.2,1
4.9,3.1,1.5,0.1,1
5.0,3.2,1.2,0.2,1
5.5,3.5,1.3,0.2,1
4.9,3.1,1.5,0.1,1
4.4,3.0,1.3,0.2,1
5.1,3.4,1.5,0.2,1
5.0,3.5,1.3,0.3,1
4.5,2.3,1.3,0.3,1
4.4,3.2,1.3,0.2,1
5.0,3.5,1.6,0.6,1
5.1,3.8,1.9,0.4,1
4.8,3.0,1.4,0.3,1
5.1,3.8,1.6,0.2,1
4.6,3.2,1.4,0.2,1
5.3,3.7,1.5,0.2,1
5.0,3.3,1.4,0.2,1
7.0,3.2,4.7,1.4,2
6.4,3.2,4.5,1.5,2
6.9,3.1,4.9,1.5,2
5.5,2.3,4.0,1.3,2
6.5,2.8,4.6,1.5,2
5.7,2.8,4.5,1.3,2
6.3,3.3,4.7,1.6,2
4.9,2.4,3.3,1.0,2
6.6,2.9,4.6,1.3,2
5.2,2.7,3.9,1.4,2
5.0,2.0,3.5,1.0,2
5.9,3.0,4.2,1.5,2
6.0,2.2,4.0,1.0,2
6.1,2.9,4.7,1.4,2
5.6,2.9,3.6,1.3,2
6.7,3.1,4.4,1.4,2
5.6,3.0,4.5,1.5,2
5.8,2.7,4.1,1.0,2
6.2,2.2,4.5,1.5,2
5.6,2.5,3.9,1.1,2
5.9,3.2,4.8,1.8,2
6.1,2.8,4.0,1.3,2
6.3,2.5,4.9,1.5,2
6.1,2.8,4.7,1.2,2
6.4,2.9,4.3,1.3,2
6.6,3.0,4.4,1.4,2
6.8,2.8,4.8,1.4,2
6.7,3.0,5.0,1.7,2
6.0,2.9,4.5,1.5,2
5.7,2.6,3.5,1.0,2
5.5,2.4,3.8,1.1,2
5.5,2.4,3.7,1.0,2
5.8,2.7,3.9,1.2,2
6.0,2.7,5.1,1.6,2
5.4,3.0,4.5,1.5,2
6.0,3.4,4.5,1.6,2
6.7,3.1,4.7,1.5,2
6.3,2.3,4.4,1.3,2
5.6,3.0,4.1,1.3,2
5.5,2.5,4.0,1.3,2
5.5,2.6,4.4,1.2,2
6.1,3.0,4.6,1.4,2
5.8,2.6,4.0,1.2,2
5.0,2.3,3.3,1.0,2
5.6,2.7,4.2,1.3,2
5.7,3.0,4.2,1.2,2
5.7,2.9,4.2,1.3,2
6.2,2.9,4.3,1.3,2
5.1,2.5,3.0,1.1,2
5.7,2.8,4.1,1.3,2
6.3,3.3,6.0,2.5,3
5.8,2.7,5.1,1.9,3
7.1,3.0,5.9,2.1,3
6.3,2.9,5.6,1.8,3
6.5,3.0,5.8,2.2,3
7.6,3.0,6.6,2.1,3
4.9,2.5,4.5,1.7,3
7.3,2.9,6.3,1.8,3
6.7,2.5,5.8,1.8,3
7.2,3.6,6.1,2.5,3
6.5,3.2,5.1,2.0,3
6.4,2.7,5.3,1.9,3
6.8,3.0,5.5,2.1,3
5.7,2.5,5.0,2.0,3
5.8,2.8,5.1,2.4,3
6.4,3.2,5.3,2.3,3
6.5,3.0,5.5,1.8,3
7.7,3.8,6.7,2.2,3
7.7,2.6,6.9,2.3,3
6.0,2.2,5.0,1.5,3
6.9,3.2,5.7,2.3,3
5.6,2.8,4.9,2.0,3
7.7,2.8,6.7,2.0,3
6.3,2.7,4.9,1.8,3
6.7,3.3,5.7,2.1,3
7.2,3.2,6.0,1.8,3
6.2,2.8,4.8,1.8,3
6.1,3.0,4.9,1.8,3
6.4,2.8,5.6,2.1,3
7.2,3.0,5.8,1.6,3
7.4,2.8,6.1,1.9,3
7.9,3.8,6.4,2.0,3
6.4,2.8,5.6,2.2,3
6.3,2.8,5.1,1.5,3
6.1,2.6,5.6,1.4,3
7.7,3.0,6.1,2.3,3
6.3,3.4,5.6,2.4,3
6.4,3.1,5.5,1.8,3
6.0,3.0,4.8,1.8,3
6.9,3.1,5.4,2.1,3
6.7,3.1,5.6,2.4,3
6.9,3.1,5.1,2.3,3
5.8,2.7,5.1,1.9,3
6.8,3.2,5.9,2.3,3
6.7,3.3,5.7,2.5,3
6.7,3.0,5.2,2.3,3
6.3,2.5,5.0,1.9,3
6.5,3.0,5.2,2.0,3
6.2,3.4,5.4,2.3,3
5.9,3.0,5.1,1.8,3];
y =  iris(:,5);
x1 = iris(:,1);
x2 = iris(:,2);
x3 = iris(:,3);
x4 = iris(:,4);
X = [x1 x2 x3 x4];




%Escolha do Número de Vizinhos;


d = 3;
K = 110;
% matriz de distancias euclidiana entre os pontos
D = euclid_dist_linalg(X,2);
%Organiza D -> para conseguir pegar construir a matriz de pesos
%que relacionará a reconstrução de cada conjunto de pontos ao
%seu vizinho.
[D_ordenado, Ind] = sort(D);


%Diferentemente do ISOMAP, nós construirremos uma matriz auxiliar
%que conterá a informação dos vizinhos, porém não com o objetivo 
%de manter a distância, mas com objetivo de saber quais serão os
%os pontos na matriz de distancias que receberão os pesos
n = size(D,1);
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

%Vamos criar a matriz de Corelação para estimação dos pesos de re
%construção do mapa, a partir dai, começa o processo de abertura
%para estimação do mapa de redimensionamento
%Minha matriz de correlação estimada será de dimensão KxK, em que
%K é a quantidade de vizinhos
alfa = 0;
beta = 0;
for i=1:n
   x_i = X(i,:)';%pegando o valor do ponto
   viz_x_i = X(Eta(i,:),:)';
   C = (viz_x_i'*viz_x_i)/(n-1);%estimando a matriz de correlação
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

%Reconstruindo o MDS clássico
M = (eye(size(W,1)) - W)'*(eye(size(W,1)) - W);
% [auto_vet, auto_val] = eigs(M, 2+1, 0);%dos menores
[Y, aut_val] = eigs(M, d+1, 0);
%autovalores, e autovetores associados, aos maiores
% [auto_vet, auto_val] = eig(M);
Y = Y(:,2:d+1)'*sqrt(n)
Y = Y';

% subplot(211)
% grid on
% hold on
% plot(Y(1:50,1),   Y(1:50,2),'b.')
% plot(Y(50:100,1), Y(50:100,2),'r.')
% plot(Y(100:150,1),Y(100:150,2),'k.')
% legend('Setosa','Versicolor','Virginica' )



subplot(212)
grid on
hold on
plot3(Y(1:50,1), Y(1:50,2), Y(1:50,3), 'b.')
plot3(Y(50:100,1), Y(50:100,2), Y(50:100,3), 'r.')
plot3(Y(100:150,1), Y(100:150,2), Y(100:150,3), 'k.')
legend('Setosa', 'Versicolor','Virginica')



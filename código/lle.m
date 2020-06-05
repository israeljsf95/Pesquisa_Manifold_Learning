%Israel Jesus Santos Filho
%LLE
clear, close all

%%Construindo o Conjunto de Pontos
n = 5;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;
mu_x = mean(X,2);
%Escolha do Número de Vizinhos;

K = 2;
% matriz de distancias euclidiana entre os pontos
D = euclid_dist_linalg(X',2);
%Organiza D -> para conseguir pegar construir a matriz de pesos
%que relacionará a reconstrução de cada conjunto de pontos ao
%seu vizinho.
[D_ordenado, Ind] = sort(D);


%Diferentemente do ISOMAP, nós construirremos uma matriz auxiliar
%que conterá a informação dos vizinhos, porém não com o objetivo 
%de manter a distância, mas com objetivo de saber quais serão os
%os pontos na matriz de distancias que receberão os pesos

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

for i=1:size(X,2)
    %Construindo a matriz de Gram
    for j = 1:n
        for k = 1:n
            G(Eta(i,j),Eta(i,k)) = (X(:,i) - X(:,Eta(i,j)))'*(X(:,i) -X(:,Eta(i,k)));
        end
    end
    G_i = inv(G+eps);
    %Usando os multiplicadores de Lagrange os pesos o ponto i
    W(i,:) = (G_i*ones(size(G,1),1))/(ones(size(G,1),1)'*(G_i*ones(size(G,1),1)));
    W(i,:) = G\ones(n,1);
end
W(isnan(W)) = 0;

%Garantindo que os pesos somem 1
for i=1:size(W,1)
    W(i,:) = W(i,:)/sum(W(i,:));
end

% 
% W(isNan(W)) = 0;







































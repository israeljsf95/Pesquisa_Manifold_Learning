%Israel ISOMAP TESTE
clear, close all
n = 1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;
ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
clf;
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled');
colormap jet(256);
%view(v1,v2);axis('equal');axis('off');
%matriz de distancias euclidianas
D = euclid_dist_linalg(X',2);
%Escolhendo os k vizinhos
k = 7;
%Organiza D 
[D_ordenado, Ind] = sort(D);
%Uma vez que D_ordenadoo está ordenado do menor para o maior, pelas colunas
%ao pegar a linha 2:k+1, pegamos os k vizinhos 
D_ordenado = D_ordenado(2:k+1,:);
Ind_ordenado = Ind(2:k+1,:);
%Este passo é extremamente importante pois já é atribuido
%através de Ind_ord e D_ord, quem são os k vizinhos de cada
%ponto presente no conjunto de dados
B = repmat(1:n, k,1);
%Convertendo os pontos no grafo de tal modo que cada elemento
%em W
W = sparse(B(:), Ind_ordenado(:), D_ordenado(:));
%Matriz que não é simétrica e as operações abaixo a transfor
%mam em uma matriz simétrica 
D_w = full(W);%D_w em que cada cada coluna contem a distancia
%do ponto da coluna com seu k ésimo vizinho, e 0 caso não seja
%vizinho
D_w = (D_w + D_w')/2;
D_w(D_w==0) = Inf;
%Isto pega adiciona em D_w o valor do ponto para ele mesmo
D_w = D_w - diag(diag(D_w));
%Algoritmo de Floyd?
for i=1:n
    Aux = repmat(D_w(:,i), 1, n)+repmat(D_w(i,:), n, 1);
    D_w = min(D_w,Aux); 
end
D_w((D_w == inf)|(isnan(D_w))) = 0;
%Reconstruindo o MDS clássico
A = (-1/2)*D_w.^2;  
H2 = eye(n) - (1/n)*(ones(n,n));
B2 = (H2*A*H2);
[auto_vet, auto_val] = eig(B2);
[aut_val_desc, ind_aut_val] = sort(auto_val, 'descend');
s = struct();
euc_Y = struct();
for i=1:10
    a_vet = auto_vet(1:end,1:i);
    a_val = auto_val(1:i,1:i);
%     a_vet = auto_vet(1:end,end:end-i+1);
%     a_val = auto_val(end:end-i+1,end:end-i+1);
    s(i).Y = a_vet*(abs(a_val).^0.5);
    euc_Y(i).Y = euclid_dist_linalg(abs(s(i).Y),1);
    J(i) = sum(sum((2*D-euc_Y(i).Y).^2))/n^2;
end
figure
plot(s(2).Y(:,1), s(2).Y(:,2),'o')
figure
scatter3(s(3).Y(:,1), s(3).Y(:,2), s(3).Y(:,3))
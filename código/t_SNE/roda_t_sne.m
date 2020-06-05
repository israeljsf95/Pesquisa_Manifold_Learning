clear all, clc, close all

%Carregando os dados
X = load('digit_xtest.csv');
Y = load('digit_ytest.csv');
X = X(1:2000,:);
Y_label = Y(1:2000);

% t = linspace(0,2*pi, 800);
% x = sin(t) + 2 * sin(2*t);
% y = cos(t) - 2*cos(2*t);
% z = -sin(3*t);
% figure
% 
% X = [x;y;z];
% X = X';

colormap jet(256)
% scatter3(x,y,z, 50, t, 'filled')

% X = load('duasDensidades.txt');
% Y_label = ones(size(X,1),1);
% Y_label(1:2000) = 1;
% Y_label(2001:end) = 2;

% Y_label = repmat([1:10], [40 1]);
% Y_label = Y_label(:);
% load('olivettifaces.mat');
% X = faces';
% %Pre-processamento: Normalização e Decomposição PCA
% Normalização
X = X - min(X(:));
X = X / max(X(:));
X = bsxfun(@minus, X, mean(X, 1));


% Cov_X = X'*X;%%matriz de covariancia de caracter?sticas
% [U2,S2,V2] = svd(X);
% n_prin_comp = 30;
% X_prin2 = U2(:,1:n_prin_comp)*S2(1:n_prin_comp,1:n_prin_comp);%%Projetando nas 30 componentes principais de X
% X_prin2 = bsxfun(@minus, X_prin2, mean(X_prin2, 1)); Não senti diferença
% ao usar esse comando 

% X_prin = bsxfun(@minus, X_prin, mean(X_prin, 1));
%PCA
% initial_dims = 30;
% if size(X, 2) < size(X, 1)  
%     C = X' * X;
% else
%     C = (1 / size(X, 1)) * (X * X');
% end
% [M, lambda] = eig(C);
% [lambda, ind] = sort(diag(lambda), 'descend');
% M = M(:,ind(1:initial_dims));
% lambda = lambda(1:initial_dims);
% if ~(size(X, 2) < size(X, 1))
%     M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
% end
% X = bsxfun(@minus, X, mean(X, 1)) * M;
% clear M lambda ind


%Calculando a matriz de Probabilidade do Kernel Gaussiano
D = pdist2(X,X).^2; %Funcionou melhor que a minha :(
% D = pdist2(X_prin2,X_prin2).^2; %Funcionou melhor que a minha :(
N = size(D,1);
P = zeros(N, N);
tol = 1e-5;
PERP = 45;
for i = 1:N
    P(i,:) = buscar_sigma(D(i,:), i, PERP, tol);
end
clear D
disp('Matriz de Probabilidades Calculada!!!');
%Garantindo que seja simétrica e que as linhas somem 1
P  = (P + P')/2;
P = max(P ./ sum(P(:)), realmin);
P = P*4;%exageração antecipada pelo algoritmo do hinton
PP = P(:).*log(P(:));
% PP(isnan(PP)) = 0;
const = sum(PP)/4;
%P = P./sum(P,2);

n_dim = 2;
%Criando os vetores de Projecao
Y = randn(N, n_dim);
% Y = 0.0001*randn(N, n_dim);
n_iter_max = 1000;
Custo = zeros(n_iter_max, 1);

%Parte Relacionada ao processo de iteração com taxa de Aprendizado
%adaptável de acordo com o artigo do Jacob, recomendado pelo hinton

eta = 500;
momentum = 0.5;
momentum_final = 0.85;
iter_troca_momentum = 0.25*n_iter_max;
parar_mod_P = 0.1*n_iter_max;
ganho_minimo = 0.01;
%Os pontos e a taxa de aprendizado devem spoder variar ao longo das
%iteracoes

incre_y = zeros(size(Y));
ganhos_eta = ones(size(Y));
iter = 1;
for i = 1:n_iter_max
    
    D_y = pdist2(Y,Y).^2;
    q = 1./(1 + D_y);
%     for i = 1:N
%         q(i,i) = 0;
%     end
%Jeito mais rápido de modificar a diagonal
    q(1:N+1:end) = 0;
    Q = max(q ./ sum(q(:)), realmin); %soma vertical para garantir que some 1
    L = (P-Q).*q;
    y_grad = 4 * (diag(sum(L,1)) - L)*Y;
    
    CC = P(:).*log(Q(:));
%     CC(isnan(CC)) = 0;
    Custo(iter) = const - sum(CC);
    
    %ganhos adaptativos
    ganhos_eta = (ganhos_eta + .05) .* (sign(y_grad) ~= sign(incre_y))...
                 + (ganhos_eta + .02) .* (sign(y_grad) == sign(incre_y)); 
    ganhos_eta(ganhos_eta < ganho_minimo) = ganho_minimo;
    incre_y = -eta*(ganhos_eta.*y_grad) + momentum*incre_y;
    Y = Y + incre_y;
%     Y = bsxfun(@minus, Y, mean(Y, 1)); 
    
%     if i == 0.35*n_iter_max
    if i == iter_troca_momentum + 1
        momentum = momentum_final;
    end
    
%     if i == 0.5*n_iter_max
    if i == parar_mod_P + 1
        P = P./4;
    end
%     fprintf('Custo iteracao %d: %.4f \n', iter, Custo(iter));
    if i
        fprintf('Custo iteracao %d: %.4f \n', iter, Custo(iter));
        scatter(Y(:,1), Y(:,2), 25, Y_label, 'filled');
%         scatter(Y(:,1), Y(:,2), 50, t, 'filled')   
        axis tight
        axis off
        drawnow
        pause(0.1)
    end
    
    iter = iter + 1;
end


scatter(Y(:,1), Y(:,2), 10, Y_label, 'filled');






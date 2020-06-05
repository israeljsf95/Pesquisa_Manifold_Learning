
clear all, clc
%Começando a rodar o programa para a MNIST


X = load('digit_xtest.csv');
X = X(1:2000,:);
Y_label = load('digit_ytest.csv');
Y_label = Y_label(1:2000);
%Como cada foto está codificado como Linha
%Pre_processando 
n_pca = 30;
X_pca = pre_proX(X, n_pca); 
fprintf('O tamanho da matriz X pre processada: \n dim_1 = %d \n dim_2 = %d \n', size(X_pca, 1), size(X_pca,2));

perp_nom = 40;
iter_max = 1000;
tol = 1e-4;
P = Matriz_prob_i(X_pca', perp_nom, tol, iter_max);
P = P*4;
disp('Matriz P Construida!!!');

Y = make_t_SNE_isr(P, 0.0002, 1e-4, 2, iter_max);
disp('Acabei!!!');

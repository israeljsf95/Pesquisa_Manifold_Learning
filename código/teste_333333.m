%%
clear all, clc, close all

n = 1000;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;
X = X/max(X(:));

ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
figure(1)
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);



figure(2)
epoch = 6500;
alfa = 1e-3;
eta = 0.3;
lambda = 0.55;
J = zeros(epoch, 8);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
W1 = -0.5 + rand(10, size(X,1));
W11 = -0.5 + rand(3, size(W1,1));
W2 = -0.5 + rand(5, size(W1,1));
W22 = -0.5 + rand(10, size(W2,1));
W3 = -0.5 + rand(3, size(W22,2));
W33 = -0.5 + rand(5, size(W3,1));
W4 = -0.5 + rand(2, size(W33,2));
W44 = -0.5 + rand(3, size(W4,1));
W5 = -0.5 + rand(3,size(W4,1));
W55 = -0.5 + rand(2, size(W5,1));
W6 = -0.5 + rand(5, size(W5,1));
W66 = -0.5 + rand(3, size(W6,1));
W7 = -0.5 + rand(10, size(W6,1));
W77 = -0.5 + rand(5, size(W7,1));
W8 = -0.5 + rand(3, size(W7,1));
W88 = -0.5 + rand(10, size(W8,1));
%Construção dos Bias
B1 = -0.5 + rand(size(W1,1),1);
B11 = -0.5 + rand(size(W11,1),1);
B2 = -0.5 + rand(size(W2,1),1);
B22 = -0.5 + rand(size(W22,1),1);
B3 = -0.5 + rand(size(W3,1),1);
B33 = -0.5 +  rand(size(W33,1),1);
B4 = -0.5 +  rand(size(W4,1),1);
B44 = -0.5 + rand(size(W44,1),1);
B5 = -0.5 + rand(size(W5,1),1);
B55 = -0.5 + rand(size(W55,1),1);
B6 = -0.5 + rand(size(W6,1),1);
B66 = -0.5 + rand(size(W66,1),1);
B7 = -0.5 + rand(size(W7,1),1);
B77 = -0.5 + rand(size(W77,1),1);
B8 = -0.5 + rand(size(W8,1),1);
B88 = -0.5 + rand(size(W88,1),1);
e = 0;
figure
subplot(811)
%Primeira camada
ee = 1;
for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE1
        a1 = tanh(W1*x + B1);
        y1 = tanh(W11*a1 + B11);

        %backprop AE1
        a1_p = 1 - a1.^2;
        y1_p = 1 - y1.^2;

        e1 = x - y1;
        delta_11 = e1.*y1_p;
        dE1_dW11 = delta_11*a1';
        dE1_dB11 = delta_11;
        delta_1  = (W11'*delta_11).*a1_p;
        dE1_dW1 = delta_1*x';
        dE1_dB1 = delta_1;
        e = e + norm(e1);
        if norm(e1) > 0.05
            W1  = W1 + alfa*dE1_dW1;
            W11 = W11 + alfa*dE1_dW11;
            B1  = B1 + alfa*dE1_dB1;
            B11 = B11 + alfa*dE1_dB11;
        else
            W1  = W1 + eta*dE1_dW1;
            W11 = W11 + eta*dE1_dW11;
            B1  = B1 + eta*dE1_dB1;
            B11 = B11 + eta*dE1_dB11;
        end
    end
    J(ee,1) = e;
    drawnow
    plot(J(:,1));
    e = 0;
    ee = ee + 1;
end
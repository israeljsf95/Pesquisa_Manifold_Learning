%%
clear all, clc, close all

n = 1500;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;
X = X/max(abs(X(:)));

ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
figure(1)
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

epoch = 3000;
alfa = 3e-3;
eta = 0.05;
lambda = 0.05;
J = zeros(epoch, 8);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
W1 = rand(15, size(X,1));
W1 = W1*sqrt(2/size(W1,1));
W11 =rand(3, size(W1,1));
W11 = W11*sqrt(2/size(W11,1));
W2 = rand(2, size(W1,1));
W2 = W2*sqrt(2/size(W2,1));
W22 =rand(15, size(W2,1));
W22 = W22*sqrt(2/size(W22,1));
% W3 = -0.5 + rand(3, size(W2,1));
% W3 = W3*sqrt(2/size(W3,1));
% W33 = -0.5 + rand(5, size(W3,1));
% W33 = W33*sqrt(2/size(W33,1));
% W4 = -0.5 + rand(2, size(W3,1));
% W4 = W4*sqrt(2/size(W4,1));
% W44 = -0.5 + rand(3, size(W4,1));
% W44 = W44*sqrt(2/size(W44,1));
% W5 = -0.5 + rand(3,size(W4,1));
% W5 = W5*sqrt(2/size(W5,1));
% W55 = -0.5 + rand(2, size(W5,1));
% W55 = W55*sqrt(2/size(W55,1));
% W6 = -0.5 + rand(5, size(W5,1));
% W6 = W6*sqrt(2/size(W6,1));
% W66 = -0.5 + rand(3, size(W6,1));
% W66 = W66*sqrt(2/size(W66,1));
% W7 = -0.5 + rand(10, size(W6,1));
% W7 = W7*sqrt(2/size(W7,1));
% W77 = -0.5 + rand(5, size(W7,1));
% W77 = W77*sqrt(2/size(W77,1));
% W8 = -0.5 + rand(3, size(W7,1));
% W8 = W8*sqrt(2/size(W8,1));
% W88 = -0.5 + rand(10, size(W8,1));
% W88 = W88*sqrt(2/size(W88,1));
%Construção dos Bias
B1 =rand(size(W1,1),1);
B1 = B1*sqrt(2/size(B1,1));
B11 = rand(size(W11,1),1);
B11 = B11*sqrt(2/size(B11,1));
B2 = rand(size(W2,1),1);
B2 = B2*sqrt(2/size(B2,1));
B22 = rand(size(W22,1),1);
B22 = B22*sqrt(2/size(B22,1));
% B3 = -0.5 + rand(size(W3,1),1);
% B3 = B3*sqrt(2/size(B3,1));
% B33 = -0.5 +  rand(size(W33,1),1);
% B33 = B33*sqrt(2/size(B33,1));
% B4 = -0.5 +  rand(size(W4,1),1);
% B4 = B4*sqrt(2/size(B4,1));
% B44 = -0.5 + rand(size(W44,1),1);
% B44 = B44*sqrt(2/size(B44,1));
% B5 = -0.5 + rand(size(W5,1),1);
% B5 = B5*sqrt(2/size(B5,1));
% B55 = -0.5 + rand(size(W55,1),1);
% B55 = B55*sqrt(2/size(B55,1));
% B6 = -0.5 + rand(size(W6,1),1);
% B6 = B6*sqrt(2/size(B6,1));
% B66 = -0.5 + rand(size(W66,1),1);
% B66 = B66*sqrt(2/size(B66,1));
% B7 = -0.5 + rand(size(W7,1),1);
% B7 = B7*sqrt(2/size(B7,1));
% B77 = -0.5 + rand(size(W77,1),1);
% B77 = B77*sqrt(2/size(B77,1));
% B8 = -0.5 + rand(size(W8,1),1);
% B8 = B8*sqrt(2/size(B8,1));
% B88 = -0.5 + rand(size(W88,1),1);
% B88 = B88*sqrt(2/size(B88,1));
e = 0;
figure
%Primeira camada
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
%         if norm(e1) > 0.05
            W1  = W1 + alfa*dE1_dW1;
            W11 = W11 + alfa*dE1_dW11;
            B1  = B1 + alfa*dE1_dB1;
            B11 = B11 + alfa*dE1_dB11;
%         else
%             W1  = W1 + eta*dE1_dW1;
%             W11 = W11 + eta*dE1_dW11;
%             B1  = B1 + eta*dE1_dB1;
%             B11 = B11 + eta*dE1_dB11;
%         end
    end
    J(ee,1) = e;
    drawnow
    plot(J(:,1));
    e = 0;
end

for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE2
        a1 = tanh(W1*x + B1);
        a2 = tanh(W2*a1 + B2);
        y2 = tanh(W22*a2 + B22);
        
        %backprop AE2
        a2_p = 1 - a2.^2;
        y2_p = 1 - y2.^2;

        e2 = a1 - y2;
        delta_22 = e2.*y2_p;
        dE1_dW22 = delta_22*a2';
        dE1_dB22 = delta_22;
        delta_2  = (W22'*delta_22).*a2_p;
        dE1_dW2 = delta_2*a1';
        dE1_dB2 = delta_2;
        e = e + norm(e2);
%         if norm(e2) > 0.05
            W2  = W2 + alfa*dE1_dW2;
            W22 = W22 + alfa*dE1_dW22;
            B2  = B2 + alfa*dE1_dB2;
            B22 = B22 + alfa*dE1_dB22;
%         else
%             W2  = W2 + eta*dE1_dW2;
%             W22 = W22 + eta*dE1_dW22;
%             B2  = B2 + eta*dE1_dB2;
%             B22 = B22 + eta*dE1_dB22;
%         end
    end
    J(ee,2) = e;
    drawnow
    plot(J(:,2));
    e = 0;
end
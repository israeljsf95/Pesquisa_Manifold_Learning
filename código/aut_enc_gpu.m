%teste com gpu

clear all, clc

n = 1500;  
% n = 5*1000; 
x = gpuArray(rand(2,n));
v = 3*pi/2*(0.1 + 2*x(1,:));
X = gpuArray(zeros(3,n));
X(1,:) = -cos(v).*v;
X(2,:) = 20*x(2,:);
X(3,:) = sin(v).*v;
X = X/norm(X);
X = [X; ones(1,size(X,2))];
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
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

%Escopo do AutoEncoder Estocástico
clear x
epoch = gpuArray(30000);
alfa = gpuArray(1e-6);
J = gpuArray(zeros(epoch, 1));
%Construção considerando os dados alocados vetores coluna
W1 = gpuArray(rand(2, size(X,1)));%aqui escolho qual dimensão os pontos serão projetados
W2 = gpuArray(rand(2, size(W1,1)));
W3 = gpuArray(rand(size(X,1),size(W2,1)));
B1 = gpuArray(rand(size(W1,1),1));
B2 = gpuArray(rand(size(W2,1),1));
B3 = gpuArray(rand(size(W3,1),1));

%Gradiente estocástico
for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,i);
        z1 = arrayfun(@tanh,(W1*x + B1));
        z2 = arrayfun(@tanh,(W2*z1 + B2));
        y  = arrayfun(@tanh,(W3*z2 + B3));
        %backpropagation
        z1_p = arrayfun(@tanh_primeg,z1);
        z2_p = arrayfun(@tanh_primeg,z2);
        y_p = array_fun(@tanh_primeg,y);
        e = y - x;
        delta_1 = e.*y_p;
        dE_dW3 = delta_1*z2';
        dE_dB3 = delta_1;
        delta_2 = (W3'*(delta_1)).*z2_p;
        dE_dW2 = delta_2*z1';
        dE_dB2 = delta_2;
        delta_3 = (W2'*(delta_2)).*z1_p;
        dE_dW1 = delta_3*x';
        dE_dB1 = delta_3;
        W1 = W1 - alfa*dE_dW1;
        B1 = B1 - alfa*dE_dB1;
        W2 = W2 - alfa*dE_dW2;
        B2 = B2 - alfa*dE_dB2;
        W3 = W3 - alfa*dE_dW3;
        B3 = B3 - alfa*dE_dB3;
    end
    J(ee) = norm(e);
    plot(J)
    drawnow
end

yy = gpuArray(zeros(size(W1,1),size(X,2)));
for i = 1:size(yy,2)
yy(:,i) = tanh(W2*tanh(W1*X(:,i) + B1) + B2);
end
figure
scatter(yy(1,:),yy(2,:))

% scatter(yy(1,1:50),yy(2,1:50), 'b')
% hold on
% scatter(yy(1,51:100),yy(2,51:100), 'r')
% scatter(yy(1,101:150),yy(2,101:150), 'g')



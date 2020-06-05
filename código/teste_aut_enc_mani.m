clear all, clc

n = 800;  
% n = 5*1000; 
x = rand(2,n);
v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
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
epoch = 30000;%epocas de treinamento
eta  = 1e-5;%taxa de aprendizado
alfa = 0.5;
W1 = rand(3, size(X,1));
W2 = rand(size(X,1), size(W1,1));
J = zeros(epoch,1);
%%Processo ded treinamento
figure
for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:size(X,2)
        x = X(:,I(i)); %amostra escolhida
        z1 = tanh(W1*x);
        %passo forward
        y = tanh(W2*z1);
        e = (x-y);
        %back_prop
        %calculo das derivadas
        z1_p = tanh_prime(z1);
        y_p  = tanh_prime(y);
        %Escrevendo o produto em função dos termos de suas derivadas!!!
        delta_1 = e.*y_p;
        dE_dW2 = delta_1*z1';
        delta_2 = (W2'*delta_1).*z1_p;
        dE_dW1 = delta_2*x';
        %atualização
        W1 = W1 + eta.*dE_dW1;
%         W2 = W2 - alfa.*dE_dW2*norm(e);
%         W3 = W3 - alfa.*dE_dW3*norm(e);
        W2 = W2 + eta.*dE_dW2;
    end
    J(ee) = norm(e);
    plot(J)
    drawnow
end



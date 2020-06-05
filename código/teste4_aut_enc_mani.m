clear all, clc, close all

n = 2000;  
m = 4;
% n = 5*1000; 
% x = rand(2,n);
x = rand(2,n)*m;
% v = 3*pi/2*(0.1 + 2*x(1,:));
X = zeros(3,n);
X(1,:) = cos(x(2,:)).*x(2,:);
X(2,:) = sin(x(2,:)).*x(2,:);
X(3,:) = x(1,:);
%pre processamento para uma rede do tipo sigmoid
X = X + abs(min(X(:)));
% X = X - repmat(mean(X,2), 1,size(X,2));
% X = X./repmat(std(X')',1,size(X,2));
X = X/max(X(:));
ms = 50;
lw = 1.5;

figure(1)
scatter3(X(1,:),X(2,:),X(3,:), ms, X(1,:), 'filled')
colormap jet(256)


epoch = 65000;


alfa1 = 1e-4;
beta1 = 2e-3;

alfa2 = 1e-5;
beta2 = 2e-4;

alfa3 = 1e-5;
beta3 = 2e-4;


alfa6 = 1e-5;
beta6 = 2e-4;

lambda = 0.0005;
J = zeros(epoch, 1);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
%AE1 parametros
W1 = 2*rand(50, size(X,1)) - 1;
W1 = W1*sqrt(2/size(W1,1));
mmt1 = zeros(size(W1));
W11 = 2*rand(50, size(W1,1)) - 1;
W11 = W11*sqrt(2/size(W11,1));
mmt11 = zeros(size(W11));

B1 = 2*rand(size(W1,1),1) - 1;
mmtB1 = zeros(size(B1));
B1 = B1*sqrt(2/size(B1,1));
B11 = 2*rand(size(W11,1),1) - 1;
mmtB11 = zeros(size(B11));
B11 = B11*sqrt(2/size(B11,1));


%AE2 parametros
W2 = 2*rand(3, size(W1,1)) - 1;
mmt2 = zeros(size(W2));
W2 = W2*sqrt(2/size(W2,1));
W22 = 2*rand(50, size(W2,1)) - 1;
W22 = W22*sqrt(2/size(W22,1));
mmt22 = zeros(size(W22));

B2 = 2*rand(size(W2,1),1) - 1;
mmtB2 = zeros(size(B2));
B2 = B2*sqrt(2/size(B2,1));
B22 = 2*rand(size(W22,1),1) - 1;
mmtB22 = zeros(size(B22));
B22 = B22*sqrt(2/size(B22,1));

%AE3 parametros

W3 = 2*rand(3, size(W22,1)) - 1;
mmt3 = zeros(size(W3));
W3 = W3*sqrt(2/size(W3,1));
B3 = 2*rand(size(W3,1),1) - 1;
mmtB3 = zeros(size(B3));
B3 = B3*sqrt(2/size(B3,1));



ee = 0;
% figure
%Primeira camada
% figure(2)
% figure(3)

%Configuração do objeto para gravar o video 


tic
for eee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE1
        a1 = relu(W1*x + B1);
        a11 = relu(W11*a1 + B11);
        a2 = relu(W2*a11 + B2);
        a22 = relu(W22*a2 + B22);
        y1 = sigmoid(W3*a22 + B3);
        
        %backprop AE1
        a1_p  = a1 > 0;
        a11_p = a11 > 0;
        a2_p  = a2 > 0;
        a22_p = a22 > 0;
        y1_p  = y1.*(1 - y1);

        e = x - y1;
        delta_1 = e.*y1_p;
        
        e1 = W3'*delta_1;
        delta_2 = a22_p.*e1;
        
        e2 = W22'*delta_2;
        delta_3 = a2_p.*e2;
        
        e3 = W2'*delta_3;
        delta_4 = a11_p.*e3;
        
        e4 = W11'*delta_4;
        delta_5 = a1_p.*e4;
        
        %atualizacao das sinapses
        dE1_dW3 = alfa1*delta_1*a22';
        dE1_dW22 = alfa1*delta_2*a2';
        dE1_dW2 = alfa1*delta_3*a11';
        dE1_dW11 = alfa1*delta_4*a1';
        dE1_dW1 = alfa1*delta_5*x';
        
        mmt3 = dE1_dW3 + beta1*mmt3;
        mmt22 = dE1_dW22 + beta1*mmt22;
        mmt2 = dE1_dW2 + beta1*mmt2;
        mmt11 = dE1_dW11 + beta1*mmt11;
        mmt1 = dE1_dW1 + beta1*mmt1;
        
        %atualizacao dos limiares
        dE1_dB3 = alfa1*delta_1;
        dE1_dB22 = alfa1*delta_2;
        dE1_dB2 = alfa1*delta_3;
        dE1_dB11 = alfa1*delta_4;
        dE1_dB1 = alfa1*delta_5;
        
        mmtB3 = dE1_dB3 + beta1*mmtB3;
        mmtB22 = dE1_dB22 + beta1*mmtB22;
        mmtB2 = dE1_dB2 + beta1*mmtB2;
        mmtB11 = dE1_dB11 + beta1*mmtB11;
        mmtB1 = dE1_dB1 + beta1*mmtB1;
        
        %visualizacao 
%         if mod(i,450) == 0
%             for j = 1:size(X,2)
%                 a111 = relu(W1*X(:,j) + B1);
%                 a222 = relu(W11*a111 + B11);
%                 a333 = relu(W2*a222 + B2);
%                 y33(:,j) = a333;
%             end
%             figure(3), subplot(211)
%             scatter3(y33(1,:),y33(2,:),y33(3,:), ms, X(1,:), 'filled')
%             colormap jet(256)
%             for j = 1:size(X,2)
%                 a111 = relu(W1*X(:,j) + B1);
%                 a222 = relu(W11*a111 + B11);
%                 a333 = relu(W2*a222 + B2);
%                 a444 = relu(W22*a333 + B22);
%                 y33(:,j) = sigmoid(W3*a444 + B3);
%             end
%             figure(3), subplot(212)
%             scatter3(y33(1,:),y33(2,:),y33(3,:), ms, X(1,:), 'filled')
%             colormap jet(256)
%             drawnow
%         end
%         
        
        
        W3  = W3 + (mmt3);
        W22  = W22 + (mmt22);
        W2  = W2 + (mmt2);
        W11 = W11 + mmt11;
        W1  = W1 + (mmt1);
        B3 = B3 + mmtB3;
        B22 = B22 + mmtB22;
        B2 = B2 + mmtB2;
        B11 = B11 + mmtB11;
        B1 = B1 + mmtB1;
        ee = ee + sum(e.^2)/3;

    end
    if mod(eee, 250) == 0
        plot(J(:,1));
    end
    J(eee,1) = ee/n;
%     drawnow    
    ee = 0;
end
toc
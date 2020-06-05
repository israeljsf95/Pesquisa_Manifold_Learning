%%
clear all, clc, close all

n = 4500;  
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


epoch = 15000;


alfa1 = 1e-5;
beta1 = 2e-4;

alfa2 = 1e-5;
beta2 = 2e-4;

alfa3 = 1e-5;
beta3 = 2e-4;


alfa6 = 1e-5;
beta6 = 2e-4;

lambda = 0.0005;
J = zeros(epoch, 2);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
%AE1 parametros
W1 = 2*rand(50, size(X,1)) - 1;
W1 = W1*sqrt(2/size(W1,1));
mmt1 = zeros(size(W1));
W11 = 2*rand(3, size(W1,1)) - 1;
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
W3 = 2*rand(50, size(W2,1)) - 1;
mmt3 = zeros(size(W3));
W3 = W3*sqrt(2/size(W3,1));
W33 = 2*rand(50, size(W3,1)) - 1;
W33 = W33*sqrt(2/size(W33,1));
mmt33 = zeros(size(W33));

B3 = 2*rand(size(W3,1),1) - 1;
mmtB3 = zeros(size(B3));
B3 = B3*sqrt(2/size(B3,1));
B33 = 2*rand(size(W33,1),1) - 1;
mmtB33 = zeros(size(B33));
B33 = B33*sqrt(2/size(B33,1));

%AE4 parametros
W4 = 2*rand(3, size(W3,1)) - 1;
mmt4 = zeros(size(W4));
W4 = W4*sqrt(2/size(W4,1));
W44 = 2*rand(50, size(W4,1)) - 1;
W44 = W44*sqrt(2/size(W44,1));
mmt44 = zeros(size(W44));

B4 = 2*rand(size(W4,1),1) - 1;
mmtB4 = zeros(size(B4));
B4 = B4*sqrt(2/size(B4,1));
B44 = 2*rand(size(W44,1),1) - 1;
mmtB44 = zeros(size(B44));
B44 = B44*sqrt(2/size(B44,1));


%AE5 parametros
W5 = 2*rand(10, size(W4,1)) - 1;
mmt5 = zeros(size(W5));
W5 = W5*sqrt(2/size(W5,1));
W55 = 2*rand(20, size(W5,1)) - 1;
W55 = W55*sqrt(2/size(W55,1));
mmt55 = zeros(size(W55));

B5 = 2*rand(size(W5,1),1) - 1;
mmtB5 = zeros(size(B5));
B5 = B5*sqrt(2/size(B5,1));
B55 = 2*rand(size(W55,1),1) - 1;
mmtB55 = zeros(size(B55));
B55 = B55*sqrt(2/size(B55,1));

%AE6 parametros
W6 = 2*rand(5, size(W5,1)) - 1;
mmt6 = zeros(size(W6));
W6 = W6*sqrt(2/size(W6,1));
W66 = 2*rand(10, size(W6,1)) - 1;
W66 = W66*sqrt(2/size(W66,1));
mmt66 = zeros(size(W66));

B6 = 2*rand(size(W6,1),1) - 1;
mmtB6 = zeros(size(B6));
B6 = B6*sqrt(2/size(B6,1));
B66 = 2*rand(size(W66,1),1) - 1;
mmtB66 = zeros(size(B66));
B66 = B66*sqrt(2/size(B66,1));

%AE7 parametros
W7 = 2*rand(3, size(W6,1)) - 1;
mmt7 = zeros(size(W7));
W7 = W7*sqrt(2/size(W7,1));
W77 = 2*rand(5, size(W7,1)) - 1;
W77 = W77*sqrt(2/size(W77,1));
mmt77 = zeros(size(W77));

B7 = 2*rand(size(W7,1),1) - 1;
mmtB7 = zeros(size(B7));
B7 = B7*sqrt(2/size(B7,1));
B77 = 2*rand(size(W77,1),1) - 1;
mmtB77 = zeros(size(B77));
B77 = B77*sqrt(2/size(B77,1));

%AE8 parametros
W8 = 2*rand(2, size(W7,1)) - 1;
mmt8 = zeros(size(W8));
W8 = W8*sqrt(2/size(W8,1));
W88 = 2*rand(3, size(W8,1)) - 1;
W88 = W88*sqrt(2/size(W88,1));
mmt88 = zeros(size(W88));

B8 = 2*rand(size(W8,1),1) - 1;
mmtB8 = zeros(size(B8));
B8 = B8*sqrt(2/size(B8,1));
B88 = 2*rand(size(W88,1),1) - 1;
mmtB88 = zeros(size(B88));
B88 = B88*sqrt(2/size(B88,1));

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
ee = 0;
% figure
%Primeira camada
subplot(811)
for eee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE1
        a1 = relu(W1*x + B1);
%         a1(a1 < 0.2) = 0;
        y1 = sigmoid(W11*a1 + B11);
%         y1(y1 < 0.3) = 0;
        
        %backprop AE1
        a1_p = a1>0;
        y1_p = y1.*(1-y1);

        e = x - y1;
        delta = e.*y1_p;
        
        e1 = W11'*delta;
        delta_1 = a1_p.*e1;
        
        %atualizacao das sinapses
        dE1_dW11 = alfa1*delta*a1';
        mmt11 = dE1_dW11 + beta1*mmt11;
        dE1_dW1 = alfa1* delta_1*x';
        mmt1 = dE1_dW1 + beta1*mmt1;
       
        %atualizacao dos limiares
        dE1_dB11 = alfa1*delta;
        mmtB11 = dE1_dB11 + beta1*mmtB11;
        dE1_dB1 = alfa1*delta_1;
        mmtB1 = dE1_dB1 + beta1*mmtB1;
      
        W1  = W1 + (mmt1);
        W11 = W11 + mmt11;
        B1 = B1 + mmtB1;
        B11 = B11 + mmtB11;
        ee = ee + sum(e.^2)/3;

    end
    J(eee,1) = ee/n;
    drawnow
    plot(J(:,1));
    ee = 0;
end

% figure
subplot(812)
for eee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE2
        a1 = relu(W1*x + B1);
        a2 = relu(W2*a1 + B2);
%         a2(a2 < 0.2) = 0;
        y2 = sigmoid(W22*a2 + B22);
        
        %backprop AE2
        a2_p = a2>0;
        y2_p = y2.*(1 - y2);

        e2 = a1 - y2;
        delta_1 = e2.*y2_p;
        
        e2_1 = W22'*delta_1;
        delta_2 = a2_p.*e2_1;
        
        %Gradiente das sinapses e dos limiares
        dE_dW22 = alfa2*delta_1*a2';
        mmt22 = dE_dW22 + beta2*mmt22;
        dE_dW2 = alfa2*delta_2*a1';
        mmt2 = dE_dW2 + beta2*mmt2;
        
        dE_dB22 = alfa2*delta_1;
        mmtB22 = dE_dB22 + beta2*mmtB22;
        dE_dB2 = alfa2*delta_2;
        mmtB2 = dE_dB2 + beta2*mmtB2;
        
        
        %atualizacao dos limiares 
        W2  = W2 + mmt2;
        W22 = W22 + mmt22;
        B2  = B2 + mmtB2;
        B22 = B22 + mmtB22;

        ee = ee + sum(e2.^2)/length(e2);
    end
    J(eee,2) = ee/n;
    drawnow
    plot(J(:,2));
    ee = 0;
end

% figure
subplot(813)
for eee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE2
        a1 = tanh(W1*x + B1);
        a2 = tanh(W2*a1 + B2);
        a3 = tanh(W3*a2 + B3);
        a3(a3 < 0.2) = 0;
        y3 = tanh(W33*a3 + B33);
        
        %backprop AE2
        a3_p = (1 - a3.^2);
        y3_p = (1 - y3.^2);

        e3 = a2 - y3;
        delta_1 = e3.*y3_p;
        
        e3_1 = W33'*delta_1;
        delta_2 = a3_p.*e3_1;
        
        %Gradiente das sinapses e dos limiares
        dE_dW33 = alfa3*delta_1*a3';
        mmt33 = dE_dW33 + beta3*mmt33;
        dE_dW3 = alfa3*delta_2*a2';
        mmt3 = dE_dW3 + beta3*mmt3;
        
        dE_dB33 = alfa3*delta_1;
        mmtB33 = dE_dB33 + beta3*mmtB33;
        dE_dB3 = alfa3*delta_2;
        mmtB3 = dE_dB3 + beta3*mmtB3;
        
        
        %atualizacao dos limiares 
        W3  = W3 + mmt3;
        W33 = W33 + mmt33;
        B3  = B3 + mmtB3;
        B33 = B33 + mmtB33;

        ee = ee + sum(e3.^2)/length(e3);
    end
    J(eee,3) = ee/n;
    drawnow
    plot(J(:,3));
    ee = 0;
end

% figure
subplot(814)
for eee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE2
        a1 = tanh(W1*x + B1);
        a2 = tanh(W2*a1 + B2);
        a3 = tanh(W3*a2 + B3);
        a4 = tanh(W4*a3 + B4);
        a4(a4 < 0.2) = 0;
        y4 = tanh(W44*a4 + B44);
        
        %backprop AE2
        a4_p = (1 - a4.^2);
        y4_p = (1 - y4.^2);

        e4 = a3 - y4;
        delta_1 = e4.*y4_p;
        
        e4_1 = W44'*delta_1;
        delta_2 = a4_p.*e4_1;
        
        %Gradiente das sinapses e dos limiares
        dE_dW44 = alfa3*delta_1*a4';
        mmt44 = dE_dW44 + beta3*mmt44;
        dE_dW4 = alfa3*delta_2*a3';
        mmt4 = dE_dW4 + beta3*mmt4;
        
        dE_dB44 = alfa3*delta_1;
        mmtB44 = dE_dB44 + beta3*mmtB44;
        dE_dB4 = alfa3*delta_2;
        mmtB4 = dE_dB4 + beta3*mmtB4;
        
        
        %atualizacao dos limiares 
        W4  = W4 + mmt4;
        W44 = W44 + mmt44;
        B4  = B4 + mmtB4;
        B44 = B44 + mmtB44;

        ee = ee + sum(e4.^2)/length(e4);
    end
    J(eee,4) = ee/n;
    drawnow
    plot(J(:,4));
    ee = 0;
end
% 
% % figure
% subplot(815)
% for eee = 1:epoch
%     I = randperm(size(X,2));
%     for i = 1:length(I)
%         x = X(:,I(i));
% 
%         %passo forward sera dividido varios autoencoders para evitar
%         %o problema do gradiente ir a zero
% 
%         %AE2
%         a1 = sigmoid(W1*x + B1);
%         a2 = sigmoid(W2*a1 + B2);
%         a3 = sigmoid(W3*a2 + B3);
%         a4 = sigmoid(W4*a3 + B4);
%         a5 = sigmoid(W5*a4 + B5);
%         a5(a5 < 0.3) = 0;
%         y5 = sigmoid(W55*a5 + B55);
%         
%         %backprop AE2
%         a5_p = a5.*(1 - a5);
%         y5_p = y5.*(1 - y5);
% 
%         e5 = a4 - y5;
%         delta_1 = e5.*y5_p;
%         
%         e5_1 = W55'*delta_1;
%         delta_2 = a5_p.*e5_1;
%         
%         %Gradiente das sinapses e dos limiares
%         dE_dW55 = alfa3*delta_1*a5';
%         mmt55 = dE_dW55 + beta3*mmt55;
%         dE_dW5 = alfa3*delta_2*a4';
%         mmt5 = dE_dW5 + beta3*mmt5;
%         
%         dE_dB55 = alfa3*delta_1;
%         mmtB55 = dE_dB55 + beta3*mmtB55;
%         dE_dB5 = alfa3*delta_2;
%         mmtB5 = dE_dB5 + beta3*mmtB5;
%         
%         
%         %atualizacao dos limiares 
%         W5  = W5 + mmt5;
%         W55 = W55 + mmt55;
%         B5  = B5 + mmtB5;
%         B55 = B55 + mmtB55;
% 
%         ee = ee + sum(e5.^2)/length(e5);
%     end
%     J(eee,5) = ee/n;
%     drawnow
%     plot(J(:,5));
%     ee = 0;
% end
% 
% % figure
% subplot(816)
% for eee = 1:epoch
%     I = randperm(size(X,2));
%     for i = 1:length(I)
%         x = X(:,I(i));
% 
%         %passo forward sera dividido varios autoencoders para evitar
%         %o problema do gradiente ir a zero
% 
%         %AE2
%         a1 = sigmoid(W1*x + B1);
%         a2 = sigmoid(W2*a1 + B2);
%         a3 = sigmoid(W3*a2 + B3);
%         a4 = sigmoid(W4*a3 + B4);
%         a5 = sigmoid(W5*a4 + B5);
%         a6 = sigmoid(W6*a5 + B6);
%         a6(a6 < 0.3) = 0;
%         y6 = sigmoid(W66*a6 + B66);
%         
%         %backprop AE2
%         a6_p = a6.*(1 - a6);
%         y6_p = y6.*(1 - y6);
% 
%         e6 = a5 - y6;
%         delta_1 = e6.*y6_p;
%         
%         e6_1 = W66'*delta_1;
%         delta_2 = a6_p.*e6_1;
%         
%         %Gradiente das sinapses e dos limiares
%         dE_dW66 = alfa3*delta_1*a6';
%         mmt66 = dE_dW66 + beta3*mmt66;
%         dE_dW6 = alfa3*delta_2*a5';
%         mmt6 = dE_dW6 + beta3*mmt6;
%         
%         dE_dB66 = alfa3*delta_1;
%         mmtB66 = dE_dB66 + beta3*mmtB66;
%         dE_dB6 = alfa3*delta_2;
%         mmtB6 = dE_dB6 + beta3*mmtB6;
%         
%         
%         %atualizacao dos limiares 
%         W6  = W6 + mmt6;
%         W66 = W66 + mmt66;
%         B6  = B6 + mmtB6;
%         B66 = B66 + mmtB66;
% 
%         ee = ee + sum(e6.^2)/length(e6);
%     end
%     J(eee,6) = ee/n;
%     drawnow
%     plot(J(:,6));
%     ee = 0;
% end
% 
% subplot(817)
% for eee = 1:epoch
%     I = randperm(size(X,2));
%     for i = 1:length(I)
%         x = X(:,I(i));
% 
%         %passo forward sera dividido varios autoencoders para evitar
%         %o problema do gradiente ir a zero
% 
%         %AE2
%         a1 = sigmoid(W1*x + B1);
%         a2 = sigmoid(W2*a1 + B2);
%         a3 = sigmoid(W3*a2 + B3);
%         a4 = sigmoid(W4*a3 + B4);
%         a5 = sigmoid(W5*a4 + B5);
%         a6 = sigmoid(W6*a5 + B6);
%         a7 = sigmoid(W7*a6 + B7);
%         a7(a7 < 0.3) = 0;
%         y7 = sigmoid(W77*a7 + B77);
%         
%         %backprop AE2
%         a7_p = a7.*(1 - a7);
%         y7_p = y7.*(1 - y7);
% 
%         e7 = a6 - y7;
%         delta_1 = e7.*y7_p;
%         
%         e7_1 = W77'*delta_1;
%         delta_2 = a7_p.*e7_1;
%         
%         %Gradiente das sinapses e dos limiares
%         dE_dW77 = alfa3*delta_1*a7';
%         mmt77 = dE_dW77 + beta3*mmt77;
%         dE_dW7 = alfa3*delta_2*a6';
%         mmt7 = dE_dW7 + beta3*mmt7;
%         
%         dE_dB77 = alfa3*delta_1;
%         mmtB77 = dE_dB77 + beta3*mmtB77;
%         dE_dB7 = alfa3*delta_2;
%         mmtB7 = dE_dB7 + beta3*mmtB7;
%         
%         
%         %atualizacao dos limiares 
%         W7  = W7 + mmt7;
%         W77 = W77 + mmt77;
%         B7  = B7 + mmtB7;
%         B77 = B77 + mmtB77;
% 
%         ee = ee + sum(e7.^2)/length(e7);
%     end
%     J(eee,7) = ee/n;
%     drawnow
%     plot(J(:,7));
%     ee = 0;
% end
% 
% subplot(818)
% for eee = 1:epoch
%     I = randperm(size(X,2));
%     for i = 1:length(I)
%         x = X(:,I(i));
% 
%         %passo forward sera dividido varios autoencoders para evitar
%         %o problema do gradiente ir a zero
% 
%         %AE2
%         a1 = sigmoid(W1*x + B1);
%         a2 = sigmoid(W2*a1 + B2);
%         a3 = sigmoid(W3*a2 + B3);
%         a4 = sigmoid(W4*a3 + B4);
%         a5 = sigmoid(W5*a4 + B5);
%         a6 = sigmoid(W6*a5 + B6);
%         a7 = sigmoid(W7*a6 + B7);
%         a8 = sigmoid(W8*a7 + B8);
%         a8(a8 < 0.3) = 0;
%         y8 = sigmoid(W88*a8 + B88);
%         
%         %backprop AE2
%         a8_p = a8.*(1 - a8);
%         y8_p = y8.*(1 - y8);
% 
%         e8 = a7 - y8;
%         delta_1 = e8.*y8_p;
%         
%         e8_1 = W88'*delta_1;
%         delta_2 = a8_p.*e8_1;
%         
%         %Gradiente das sinapses e dos limiares
%         dE_dW88 = alfa6*delta_1*a8';
%         mmt88 = dE_dW88 + beta6*mmt88;
%         dE_dW8 = alfa6*delta_2*a7';
%         mmt8 = dE_dW8 + beta6*mmt8;
%         
%         dE_dB88 = alfa6*delta_1;
%         mmtB88 = dE_dB88 + beta6*mmtB88;
%         dE_dB8 = alfa6*delta_2;
%         mmtB8 = dE_dB8 + beta6*mmtB8;
%         
%         
%         %atualizacao dos limiares 
%         W8  = W8 + mmt8;
%         W88 = W88 + mmt88;
%         B8  = B8 + mmtB8;
%         B88 = B88 + mmtB88;
% 
%         ee = ee + sum(e8.^2)/length(e8);
%     end
%     J(eee,8) = ee/n;
%     drawnow
%     plot(J(:,8));
%     ee = 0;
% end
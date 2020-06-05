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
X = X/max(abs(X(:)))/10;

ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;
figure(1)
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

epoch = 1500;
alfa = 3e-3;
eta = 0.05;
lambda = 0.05;
J = zeros(epoch, 8);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
W1 = -0.5 + rand(10, size(X,1));
W1 = W1*sqrt(2/size(W1,1));
W11 = -0.5 + rand(3, size(W1,1));
W11 = W11*sqrt(2/size(W11,1));
W2 = -0.5 + rand(5, size(W1,1));
W2 = W2*sqrt(2/size(W2,1));
W22 = -0.5 + rand(10, size(W2,1));
W22 = W22*sqrt(2/size(W22,1));
W3 = -0.5 + rand(3, size(W2,1));
W3 = W3*sqrt(2/size(W3,1));
W33 = -0.5 + rand(5, size(W3,1));
W33 = W33*sqrt(2/size(W33,1));
W4 = -0.5 + rand(2, size(W3,1));
W4 = W4*sqrt(2/size(W4,1));
W44 = -0.5 + rand(3, size(W4,1));
W44 = W44*sqrt(2/size(W44,1));
W5 = -0.5 + rand(3,size(W4,1));
W5 = W5*sqrt(2/size(W5,1));
W55 = -0.5 + rand(2, size(W5,1));
W55 = W55*sqrt(2/size(W55,1));
W6 = -0.5 + rand(5, size(W5,1));
W6 = W6*sqrt(2/size(W6,1));
W66 = -0.5 + rand(3, size(W6,1));
W66 = W66*sqrt(2/size(W66,1));
W7 = -0.5 + rand(10, size(W6,1));
W7 = W7*sqrt(2/size(W7,1));
W77 = -0.5 + rand(5, size(W7,1));
W77 = W77*sqrt(2/size(W77,1));
W8 = -0.5 + rand(3, size(W7,1));
W8 = W8*sqrt(2/size(W8,1));
W88 = -0.5 + rand(10, size(W8,1));
W88 = W88*sqrt(2/size(W88,1));
%Construção dos Bias
B1 = -0.5 + rand(size(W1,1),1);
B1 = B1*sqrt(2/size(B1,1));
B11 = -0.5 + rand(size(W11,1),1);
B11 = B11*sqrt(2/size(B11,1));
B2 = -0.5 + rand(size(W2,1),1);
B2 = B2*sqrt(2/size(B2,1));
B22 = -0.5 + rand(size(W22,1),1);
B22 = B22*sqrt(2/size(B22,1));
B3 = -0.5 + rand(size(W3,1),1);
B3 = B3*sqrt(2/size(B3,1));
B33 = -0.5 +  rand(size(W33,1),1);
B33 = B33*sqrt(2/size(B33,1));
B4 = -0.5 +  rand(size(W4,1),1);
B4 = B4*sqrt(2/size(B4,1));
B44 = -0.5 + rand(size(W44,1),1);
B44 = B44*sqrt(2/size(B44,1));
B5 = -0.5 + rand(size(W5,1),1);
B5 = B5*sqrt(2/size(B5,1));
B55 = -0.5 + rand(size(W55,1),1);
B55 = B55*sqrt(2/size(B55,1));
B6 = -0.5 + rand(size(W6,1),1);
B6 = B6*sqrt(2/size(B6,1));
B66 = -0.5 + rand(size(W66,1),1);
B66 = B66*sqrt(2/size(B66,1));
B7 = -0.5 + rand(size(W7,1),1);
B7 = B7*sqrt(2/size(B7,1));
B77 = -0.5 + rand(size(W77,1),1);
B77 = B77*sqrt(2/size(B77,1));
B8 = -0.5 + rand(size(W8,1),1);
B8 = B8*sqrt(2/size(B8,1));
B88 = -0.5 + rand(size(W88,1),1);
B88 = B88*sqrt(2/size(B88,1));
e = 0;
figure
subplot(811)
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

%Segunda Camada
subplot(812)
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

subplot(813)
for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));

        %passo forward sera dividido varios autoencoders para evitar
        %o problema do gradiente ir a zero

        %AE2
        a1 = tanh(W1*x + B1);
        a2 = tanh(W2*a1 + B2);
        a3 = tanh(W3*a2 + B3);
        y3 = tanh(W33*a3 + B33);

        
        %backprop AE2
        a3_p = 1 - a3.^2;
        y3_p = 1 - y3.^2;

        e3 = a2 - y3;
        delta_33 = e3.*y3_p;
        dE1_dW33 = delta_33*a3';
        dE1_dB33 = delta_33;
        delta_3  = (W33'*delta_33).*a3_p;
        dE1_dW3 = delta_3*a2';
        dE1_dB3 = delta_3;
        e = e + norm(e3);
%         if norm(e3) > 0.05
            W3  = W3 + alfa*dE1_dW3;
            W33 = W33 + alfa*dE1_dW33;
            B3  = B3 + alfa*dE1_dB3;
            B33 = B33 + alfa*dE1_dB33;
%         else
%             W3  = W3 + eta*dE1_dW3;
%             W33 = W33 + eta*dE1_dW33;
%             B3  = B3 + eta*dE1_dB3;
%             B33 = B33 + eta*dE1_dB33;
%         end
    end
    J(ee,3) = e;
    drawnow
    plot(J(:,3));
    e = 0;
end

subplot(814)
for ee = 1:epoch
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
        y4 = tanh(W44*a4 + B44);
        
        %backprop AE2
        a4_p = 1 - a4.^2;
        y4_p = 1 - y4.^2;

        e4 = a3 - y4;
        delta_44 = e4.*y4_p;
        dE1_dW44 = delta_44*a4';
        dE1_dB44 = delta_44;
        delta_4  = (W44'*delta_44).*a4_p;
        dE1_dW4 = delta_4*a3';
        dE1_dB4 = delta_4;
        e = e + norm(e4);
%         if norm(e4) > 0.2
            W4  = W4 + 1e-3*dE1_dW4;
            W44 = W44 + 1e-3*dE1_dW44;
            B4  = B4 + 1e-3*dE1_dB4;
            B44 = B44 + 1e-3*dE1_dB44;
%         else
%             W4  = W4 + alfa*dE1_dW4;
%             W44 = W44 + alfa*dE1_dW44;
%             B4  = B4 + alfa*dE1_dB4;
%             B44 = B44 + alfa*dE1_dB44;
%         end
    end
    J(ee,4) = e;
    drawnow
    plot(J(:,4));
    e = 0;
end

subplot(815)
for ee = 1:epoch
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
        a5 = tanh(W5*a4 + B5);
        y5 = tanh(W55*a5 + B55);

        %backprop AE5
        a5_p = 1 - a5.^2;
        y5_p = 1 - y5.^2;

        e5 = a4 - y5;
        delta_55 = e5.*y5_p;
        dE1_dW55 = delta_55*a5';
        dE1_dB55 = delta_55;
        delta_5  = (W55'*delta_55).*a5_p;
        dE1_dW5 = delta_5*a4';
        dE1_dB5 = delta_5;
        e = e + norm(e5);
%         if norm(e5) > 0.05
            W5  = W5 + alfa*dE1_dW5;
            W55 = W55 + alfa*dE1_dW55;
            B5  = B5 + alfa*dE1_dB5;
            B55 = B55 + alfa*dE1_dB55;
%         else
%             W5  = W5 + eta*dE1_dW5;
%             W55 = W55 + eta*dE1_dW55;
%             B5  = B5 + eta*dE1_dB5;
%             B55 = B55 + eta*dE1_dB55;
%         end
    end
    J(ee,5) = e;
    drawnow
    plot(J(:,5));
    e = 0;
end

subplot(816)
for ee = 1:epoch
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
        a5 = tanh(W5*a4 + B5);
        a6 = tanh(W6*a5 + B6);
        y6 = tanh(W66*a6 + B66);

        %backprop AE6
        a6_p = 1 - a6.^2;
        y6_p = 1 - y6.^2;

        e6 = a5 - y6;
        delta_66 = e6.*y6_p;
        dE1_dW66 = delta_66*a6';
        dE1_dB66 = delta_66;
        delta_6  = (W66'*delta_66).*a6_p;
        dE1_dW6 = delta_6*a5';
        dE1_dB6 = delta_6;
        e = e + norm(e6);
%         if norm(e6) > 0.05
            W6  = W6 + alfa*dE1_dW6;
            W66 = W66 + alfa*dE1_dW66;
            B6  = B6 + alfa*dE1_dB6;
            B66 = B66 + alfa*dE1_dB66;
%         else
%             W6  = W6 + eta*dE1_dW6;
%             W66 = W66 + eta*dE1_dW66;
%             B6  = B6 + eta*dE1_dB6;
%             B66 = B66 + eta*dE1_dB66;
%         end
    end
    J(ee,6) = e;
    drawnow
    plot(J(:,6));
    e = 0;
end

subplot(817)
for ee = 1:epoch
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
        a5 = tanh(W5*a4 + B5);
        a6 = tanh(W6*a5 + B6);
        a7 = tanh(W7*a6 + B7);
        y7 = tanh(W77*a7 + B77);

        %backprop AE7
        a7_p = 1 - a7.^2;
        y7_p = 1 - y7.^2;
       
        e7 = a6 - y7;
        delta_77 = e7.*y7_p;
        dE1_dW77 = delta_77*a7';
        dE1_dB77 = delta_77;
        delta_7  = (W77'*delta_77).*a7_p;
        dE1_dW7 = delta_7*a6';
        dE1_dB7 = delta_7;
        e = e + norm(e7);
%         if norm(e7) > 0.05
            W7  = W7 + alfa*dE1_dW7;
            W77 = W77 + alfa*dE1_dW77;
            B7  = B7 + alfa*dE1_dB7;
            B77 = B77 + alfa*dE1_dB77;
%         else
%             W7  = W7 + eta*dE1_dW7;
%             W77 = W77 + eta*dE1_dW77;
%             B7  = B7 + eta*dE1_dB7;
%             B77 = B77 + eta*dE1_dB77;
%         end
    end
    J(ee,7) = e;
    drawnow
    plot(J(:,7));
    e = 0;
end

subplot(818)
for ee = 1:epoch
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
        a5 = tanh(W5*a4 + B5);
        a6 = tanh(W6*a5 + B6);
        a7 = tanh(W7*a6 + B7);
        a8 = tanh(W8*a7 + B8);
        y8 = tanh(W88*a8 + B88);

        %backprop AE8
        a8_p = 1 - a8.^2;
        y8_p = 1 - y8.^2;
       
        e8 = a7 - y8;
        delta_88 = e8.*y8_p;
        dE1_dW88 = delta_88*a8';
        dE1_dB88 = delta_88;
        delta_8  = (W88'*delta_88).*a8_p;
        dE1_dW8 = delta_8*a7';
        dE1_dB8 = delta_8;
        e = e + norm(e8);
%         if norm(e8) > 0.05
            W8  = W8 + alfa*dE1_dW8;
            W88 = W88 + alfa*dE1_dW88;
            B8  = B8 + alfa*dE1_dB8;
            B88 = B88 + alfa*dE1_dB88;
%         else
%             W8  = W8 + eta*dE1_dW8;
%             W88 = W88 + eta*dE1_dW88;
%             B8  = B8 + eta*dE1_dB8;
%             B88 = B88 + eta*dE1_dB88;
%         end
    end
    J(ee,8) = e;
    drawnow
    plot(J(:,8));
    e = 0;
end

%%

y22 = zeros(2, size(X,2));
y33 = zeros(3, size(X,2));

for i = 1:size(X,2)
    a1 = tanh(W1*X(:,i) + B1);
    a2 = tanh(W2*a1 + B2);
    a3 = tanh(W3*a2 + B3);
    a4 = tanh(W4*a3 + B4);
    y22(:,i) = a4;
    y33(:,i) = a3;
end




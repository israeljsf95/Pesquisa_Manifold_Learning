clear all, clc

n = 600;  
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
epoch = 2500;
alfa = 1e-4;
eta = 0.3;
lambda = 0.55;
J = zeros(epoch, 8);
%Construção considerando os dados alocados vetores coluna
%Construção dos Pesos Sinápticos
W1 = rand(5, size(X,1));
W11 = rand(3, size(W1,1));
W2 = rand(4, size(W1,1));
W22 = rand(5, size(W2,1));
W3 = rand(3, size(W22,2));
W33 = rand(4, size(W3,1));
W4 = rand(2, size(W33,2));
W44 = rand(3, size(W4,1));
W5 = rand(3,size(W4,1));
W55 = rand(2, size(W5,1));
W6 = rand(4, size(W5,1));
W66 = rand(3, size(W6,1));
W7 = rand(5, size(W6,1));
W77 = rand(4, size(W7,1));
W8 = rand(3, size(W7,1));
W88 = rand(5, size(W8,1));
%Construção dos Bias
B1 = rand(size(W1,1),1);
B11 = rand(size(W11,1),1);
B2 = rand(size(W2,1),1);
B22 = rand(size(W22,1),1);
B3 = rand(size(W3,1),1);
B33 = rand(size(W33,1),1);
B4 = rand(size(W4,1),1);
B44 = rand(size(W44,1),1);
B5 = rand(size(W5,1),1);
B55 = rand(size(W55,1),1);
B6 = rand(size(W6,1),1);
B66 = rand(size(W66,1),1);
B7 = rand(size(W7,1),1);
B77 = rand(size(W77,1),1);
B8 = rand(size(W8,1),1);
B88 = rand(size(W88,1),1);
e = 0;

for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:size(X,2)
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
        
        %AE2
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
        
        %AE3
        a3 = tanh(W3*a2 + B3);
        y3 = tanh(W33*a3 + B33);
        
        %backprop AE3
        a3_p = 1 - a3.^2;
        y3_p = 1 - y3.^2;
        
        e3 = a2 - y3;
        delta_33 = e3.*y3_p;
        dE1_dW33 = delta_33*a3';
        dE1_dB33 = delta_33;
        delta_3  = (W33'*delta_33).*a3_p;
        dE1_dW3 = delta_3*a2';
        dE1_dB3 = delta_3;
        
        %AE4
        a4 = tanh(W4*a3 + B4);
        y4 = tanh(W44*a4 + B44);
        
        %backprop AE4
        a4_p = 1 - a4.^2;
        y4_p = 1 - y4.^2;
        
        e4 = a3 - y4;
        delta_44 = e4.*y4_p;
        dE1_dW44 = delta_44*a4';
        dE1_dB44 = delta_44;
        delta_4  = (W44'*delta_44).*a4_p;
        dE1_dW4 = delta_4*a3';
        dE1_dB4 = delta_4;
        
        %AE5
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
        
        %AE6
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
        
        %AE7
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

        %AE8
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
        %erro

        W1  = W1 + alfa*dE1_dW1;
        W11 = W11 + alfa*dE1_dW11;
        W2  = W2 + alfa*dE1_dW2;
        W22 = W22 + alfa*dE1_dW22;
        W3  = W3 + alfa*dE1_dW3;
        W33 = W33 + alfa*dE1_dW33;
        W4  = W4 + alfa*dE1_dW4;
        W44 = W44 + alfa*dE1_dW44;
        W5  = W5 + alfa*dE1_dW5;
        W55 = W55 + alfa*dE1_dW55;
        W6  = W6 + alfa*dE1_dW6;
        W66 = W66 + alfa*dE1_dW66;
        W7  = W7 + alfa*dE1_dW7;
        W77 = W77 + alfa*dE1_dW77;
        W8  = W8 + alfa*dE1_dW8;
        W88 = W88 + alfa*dE1_dW88;
        %Construção dos Bias
        B1  = B1 + alfa*dE1_dB1;
        B11 = B11 + alfa*dE1_dB11;
        B2  = B2 + alfa*dE1_dB2;
        B22 = B22 + alfa*dE1_dB22;
        B3  = B3 + alfa*dE1_dB3;
        B33 = B33 + alfa*dE1_dB33;
        B4  = B4 + alfa*dE1_dB4;
        B44 = B44 + alfa*dE1_dB44;
        B5  = B5 + alfa*dE1_dB5;
        B55 = B55 + alfa*dE1_dB55;
        B6  = B6 + alfa*dE1_dB6;
        B66 = B66 + alfa*dE1_dB66;
        B7  = B7 + alfa*dE1_dB7;
        B77 = B77 + alfa*dE1_dB77;
        B8  = B8 + alfa*dE1_dB8;
        B88 = B88 + alfa*dE1_dB88;
        e = e + norm(x-a8);      
    end
    J(ee) = e;
    plot(J)
    drawnow
    e = 0;
end


yy2 = zeros(size(W4,1),size(X,2));
for i = 1:size(yy2,2)
yy2(:,i) = tanh(W4*(W3*tanh(W2*tanh(W1*X(:,i) + B1) + B2) + B3) + B4);
end

yy3 = zeros(size(W3,1),size(X,2));
for i = 1:size(yy3,2)
yy3(:,i) = tanh(W3*tanh(W2*tanh(W1*X(:,i) + B1) + B2) + B3);
end

y = zeros(size(X));
for i = 1:size(X,2)
    y(i) = tanh(W8*tanh(W7*tanh(W6*tanh(W5*tanh(W4*(W3*tanh(W2*tanh(W1*X(:,i) + B1) + B2) + B3) + B4)+B5) + B6)+ B7) + B8);
end
clear all, clc, close all

n = 1000;  
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
figure
scatter3(X(1,:),X(2,:),X(3,:), ms, v, 'filled')
colormap jet(256)
view(v1,v2);

%Escopo do AutoEncoder Estocástico
clear x
epoch = 4000;
alfa = 1e-3;
lambda = 3e-1;
J = zeros(epoch, 1);
%Construção considerando os dados alocados vetores coluna
W1 = -0.3 + rand(5, size(X,1));%aqui escolho qual dimensão os pontos serão projetados
W1 = W1*sqrt(2/size(W1,1));
W2 = -0.3 + rand(2, size(W1,1));
W2 = W2*sqrt(2/size(W2,1));
W3 = -0.3 + rand(5,size(W2,1));
W3 = W3*sqrt(2/size(W3,1));
W4 = -0.3 + rand(3,size(W3,1));
W4 = W4*sqrt(2/size(W4,1));
B1 = -0.3 + rand(size(W1,1),1);
B1 = B1*sqrt(2/size(B1,1));
B2 = -0.3 + rand(size(W2,1),1);
B2 = B2*sqrt(2/size(B2,1));
B3 = -0.3 + rand(size(W3,1),1);
B3 = B3*sqrt(2/size(B3,1));
B4 = -0.3 + rand(size(W4,1),1);
B4 = B4*sqrt(2/size(B4,1));
%Gradiente estocástico
% Norm = zeros(epoch, 8);
e = 0;
for ee = 1:epoch
    I = randperm(size(X,2));
    for i = 1:length(I)
        x = X(:,I(i));
        z1 = tanh(W1*x + B1);
        z2 = tanh(W2*z1 + B2);
        z3 = tanh(W3*z2 + B3);
        y  = tanh(W4*z3 + B4);
        %backpropagation
        z1_p = 1 - z1.^2;
        z2_p = 1 - z2.^2;
        z3_p = 1 - z3.^2;
        y_p  = 1 - y.^2;
        
        e1 = x - y;
        delta_4 = e1.*y_p*norm(e1);
        dE_dW4 = delta_4*z3';
%         dE_dW4 = dE_dW4 + (0.001)*norm(dE_dW4);
        dE_dB4 = delta_4;
        delta_3 = (W4'*delta_4).*z3_p;
        dE_dW3 = delta_3*z2';
%         dE_dW3 = dE_dW3 + (0.001)*norm(dE_dW3);
        dE_dB3 = delta_3;
        delta_2 = (W3'*(delta_3)).*z2_p;
        dE_dW2 = delta_2*z1';
%         dE_dW2 = dE_dW2 + (0.001)*norm(dE_dW2);
        dE_dB2 = delta_2;
        delta_1 = (W2'*(delta_2)).*z1_p;
        dE_dW1 = delta_1*x';
%         dE_dW1 = dE_dW1 + (0.001)*norm(dE_dW1);
        dE_dB1 = delta_1;
%         NORM = [norm(dE_dW1) norm(dE_dW2) norm(dE_dW3) norm(dE_dW4)...
%                 norm(dE_dB1) norm(dE_dB2) norm(dE_dB3) norm(dE_dB4)];
%         if ee >=1000
%             dE_dW1 = sign(dE_dW1)*norm(Norm(ee-1,:));
%             dE_dB1 = sign(dE_dB1)*norm(Norm(ee-1,:));   
%             dE_dW2 = sign(dE_dW2)*norm(Norm(ee-1,:));
%             dE_dB2 = sign(dE_dB2)*norm(Norm(ee-1,:));
%             dE_dW3 = sign(dE_dW3)*norm(Norm(ee-1,:));
%             dE_dB3 = sign(dE_dB3)*norm(Norm(ee-1,:));
%             dE_dW4 = sign(dE_dW4)*norm(Norm(ee-1,:));
%             dE_dB4 = sign(dE_dB4)*norm(Norm(ee-1,:));
%             W1 = W1 + dE_dW1;  %- (lambda/2)*norm(W1);
%             B1 = B1 + dE_dB1;  %- (lambda/2)*norm(B1);
%             W2 = W2 + dE_dW2; %- (lambda/2)*norm(W2);
%             B2 = B2 + dE_dB2; %- (lambda/2)*norm(B2);
%             W3 = W3 + dE_dW3; %- (lambda/2)*norm(W3);
%             B3 = B3 + dE_dB3; %- (lambda/2)*norm(B3);
%             W4 = W4 + dE_dW4; %- (lambda/2)*norm(W4);
%             B4 = B4 + dE_dB4;
%         else
%             W1 = W1 + alfa*dE_dW1;  %- (lambda/2)*norm(W1);
%             B1 = B1 + alfa*dE_dB1;  %- (lambda/2)*norm(B1);
%             W2 = W2 + alfa*dE_dW2; %- (lambda/2)*norm(W2);
%             B2 = B2 + alfa*dE_dB2; %- (lambda/2)*norm(B2);
%             W3 = W3 + alfa*dE_dW3; %- (lambda/2)*norm(W3);
%             B3 = B3 + alfa*dE_dB3; %- (lambda/2)*norm(B3);
%             W4 = W4 + alfa*dE_dW4; %- (lambda/2)*norm(W4);
%             B4 = B4 + alfa*dE_dB4; %- (lambda/2)*norm(B4);
%         end
        W1 = W1 + alfa*dE_dW1;  %- (lambda/2)*norm(W1);
        B1 = B1 + alfa*dE_dB1;  %- (lambda/2)*norm(B1);
        W2 = W2 + alfa*dE_dW2; %- (lambda/2)*norm(W2);
        B2 = B2 + alfa*dE_dB2; %- (lambda/2)*norm(B2);
        W3 = W3 + alfa*dE_dW3; %- (lambda/2)*norm(W3);
        B3 = B3 + alfa*dE_dB3; %- (lambda/2)*norm(B3);
        W4 = W4 + alfa*dE_dW4; %- (lambda/2)*norm(W4);
        B4 = B4 + alfa*dE_dB4; %- (lambda/2)*norm(B4);
        e = e + norm(e1);
    end
    J(ee) = norm(e);
    plot(J)
    drawnow
    e = 0;
%     Norm(ee,:) = NORM;
end

% yy = zeros(size(W1,1),size(X,2));
% for i = 1:size(yy,2)
% yy(:,i) = tanh(W2*tanh(W1*X(:,i) + B1) + B2);
% end
% figure
% scatter(yy(1,:),yy(2,:))

% scatter(yy(1,1:50),yy(2,1:50), 'b')
% hold on
% scatter(yy(1,51:100),yy(2,51:100), 'r')
% scatter(yy(1,101:150),yy(2,101:150), 'g')



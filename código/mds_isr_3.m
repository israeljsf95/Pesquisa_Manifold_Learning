function [E_1, E_2, n_dim, P2, P] = mds_isr_3(dr, iter, passo)
%dr e uma matriz simetrica de dissimiliradade e quadrada
n_dim = size(dr,1);%recebe  o n�mero de pontos, linhas da matriz dr
E_1= zeros(n_dim, iter); %% E_1 � um vetor linha que vai receber os valores da quantidade de stress
vetor_aux = 1:n_dim;%vetor_aux ajuda nas contas da derivada
repet = 10;
E_2 = zeros(n_dim,repet);
for kk = 1:repet    
    for i = 1:n_dim
        P = 100*rand(n_dim, i);%matriz que variara as dimensoes ate n. � rand para evitar a proximidade a origem.
        for j = 1:iter
            N = randi(n_dim,1,1);%sorteio de um ponto da matr   iz P;
            dp = grad_aux(P,i,N);%dp � o vetor de distancias que sera usada nas aproxima��es das dimens�es
            dn = sqrt(sum(dp.^2,1));%dn � o vetor de distancia euclidiana entre o ponto escolhido
            v = (dn - dr(N,vetor_aux~=N)).^2;%v � o vetor de erro do ponto 
            J = norm(v);%J e o estresse associado ao ponto N
            dJ_n = zeros(i,n_dim-1);%a medida em que a dimens�o cresce, dJ_n se tornara uma matriz de atualizacao dos pontos de acordo com a dimens�o do atual, indice i
            %calculo da derivada para para atualiza��o da coordenada
            for k=1:i 
                dJ_n(k,:) = (-1)*(1/J)*(((dn - dr(N,vetor_aux~=N)).*dp(k,:))./dn);
            end
            dJ_n = sum(dJ_n, 2);
            P(N,:) = P(N,:) - passo*(dJ_n)';%atualizaca��o do Ponto N
            E_1(i,j) = J;
            if i==2
                P2 = P;
            end;
        end
    end
    E = mean(E_1, 2);%tiro a media do erro por dimensao e normalizo pelo primeiro valor para fazer a curva do joelho entre 0-1 
    E  = E./E(1);
    E_2(:,kk) = E;    
end
E_2 = mean(E_1, 2);%tiro a media do erro final por dimensao e normalizo pelo primeiro valor para fazer a curva do joelho entre 0-1 
E_2  = E_2./E_2(1);

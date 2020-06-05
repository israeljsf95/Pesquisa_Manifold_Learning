%forma de calcular a matriz de distancias

function [m_dist_1] = euclid_dist_isr(P, mode)
%matriz m x n, como vou usar no mds, a quantidade de linhas e igual a
%quantidade de dados provenientes na matriz
[m, n] = size(P); 
m_dist_1 = zeros(m, m);
%m_dist_2 = zeros(m, m);
if (mode == 1)    
    for i=1:m
        for j=1:m
            m_dist_1(i,j) = sqrt(sum((P(i,:) - P(j,:)).^2));
        end
    end
elseif (mode == 2)
    for i=1:m
        for j=1:m
            m_dist_1(i,j) = sum((P(i,:) - P(j,:)));
        end
    end
else
    fprintf('erro\n');
end
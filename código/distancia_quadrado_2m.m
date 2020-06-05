function [S] = distancia_quadrado_2m(Y, X)
    %funcao que calcula a distancia euclidiana ao quadrado
    %supondo que os dados são vetores coluna em X e em Y
    S = repmat(sum(X.*X),size(Y,2),1) + repmat(sum(Y.*Y)', 1, size(X,2)) - 2*Y'*X;
end
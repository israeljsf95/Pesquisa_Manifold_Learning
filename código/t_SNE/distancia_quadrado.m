function [S] = distancia_quadrado(X)
    %funcao que calcula a distancia euclidiana ao quadrado
    %supondo que os dados são vetores coluna em X
    D2 = repmat(sum(X.*X),size(X,2),1);
    S = D2 + D2' - 2*X'*X;
end
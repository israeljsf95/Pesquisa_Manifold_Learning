
function M = euclid_dist_linalg(p, m)
%matriz m x n, como vou usar no mds, a quantidade de linhas e igual a
%quantidade de dados provenientes na matriz
    B = p*p';
    n = size(B, 1);
    %pegando os elementos diagonais de A
    c = diag(B);
    oness = ones(n, 1);
    %calculando a matriz de distancias ao quadrado
    M = c*oness' + oness*c' - 2*B;
    %tirando a raiz de P
    if m==1
        M;
    elseif m==2
        M = sqrt(M);
    else
        fprintf('erro \n');
    end
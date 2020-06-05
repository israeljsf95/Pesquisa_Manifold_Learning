function Q = criar_matriz_Q(P, Y)
    %criando a matriz de distanicas
    D_y = pdist2(Y,Y).^2;
    q = 1./(1 + D_y);
    %zerando a diagonal de q
    for i = 1:size(q,1)
        q(i,i) = 0;
    end
    %Garandtindo que as linhas de Q somem 1
    Q = q./sum(q,2);
end
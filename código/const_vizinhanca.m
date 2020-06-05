function [ND2, A] = const_vizinhanca(X, options)
% Construct neighborhood system on data.
% [ND2, A] = DATA NEIGHBORHOOD(X, OPTIONS)
% INPUTS:
% X: (D x n) matrix,
% in which X(:,i) is the ith point in R^D.
% OPTIONS:
% .epsilon(double) for epsilon-neighbors.
% .k (integer or n-vector) for k-neighbors.
% .symmetric (Boolean)
% if =1, output a symmetric pre-weight matrix.
% Otherwise, it may be asymmetric
% .verb (Boolean) if =1, display the comments.
% OUTPUTS:
% ND2 (double,sparse): sparse matrix, in which
% each row contains only the square distances
% from the i-th point to its neighbors. i.e.,
% ND2(i,j)=0, if point j is not a neighbor of
% the i-th point.
% Otherwise, ND2(i,j) = S(i,j).
% A (Boolean, sparse): the adjacency matrix.
% A(i,j) = 1, if x j is a neighbor of x i.
% Otherwise, A(i,j) = 0.
% Example.
% X = rand(200,400);
% options.k = 10;
% options.symmetric = 0;
% options.verb = 1;
% [ND2,A]= data neighborhood(X, options);
% Initialize options.
    options.null=0;
    n = size(X,2);
    emode = false;
    if isfield(options, 'verb')
        verb = options.verb;
    else
        verb = 1;
    end

    if isfield(options, 'epsilon')
        epsilon = options.epsilon;
        emode = true;
        elseif isfield(options, 'k')
            k = options.k;
    % Check the validity of k
            if (numel(k)~=1)&&(numel(k)~=n)
                error('k must be integer or n-vector');
            elseif numel(k)==1
                k=k*ones(1,n);
            else
                k=k(:);
            end
    else
        k = 10*ones(1,n);
    end
    
    if isfield(options, 'symmetric')
        symmetric = options.symmetric;
    else
        symmetric = 0;
    end
    % Compute square distance matrix.
    if verb
        disp('- Calculando a matriz de distancias.');
    end
    S = distancia_quadrado(X);
    clear X
    % Construct neighborhood and compute
    % square distances for neighbors.
    if emode % Vizinhanca bola epsillon.
        if verb
            disp('- Calculando a vizinhança na bola eps');
        end
        S(S> epsilon^2) = 0;
        ND2 = sparse(S);
    else % Calculo dos K vizinhos.
        if verb
            disp('- Calculando os K vizinhos.');
        end
        ND2 = sparse(n,n);
        [sortS, indexS] = sort(S);
        for i = 1:n
            nbri = indexS(2:k(i)+1,i);
            ND2(nbri,i) = sortS(2:k(i)+1,i);
        end
    end
    if symmetric
        ND2 = (ND2 + ND2')/2;
    end
    
    A = (ND2>0);

end
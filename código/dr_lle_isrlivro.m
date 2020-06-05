function [Y, ev] = dr_lle_isrlivro(X,d,options)
% Locally Linear Embedding for DR
% SYNTAX:
% Y = drlle(X,d,options)
% INPUTS:
% X: D x n data matrix, X(:,i) is ith point.
% d: Dimension of output data.
% OPTIONS:
% .epsilon: epsilon-neighbors.
% .k: k-neighbors.
% .nonnegativeW: True for nonnegative weights.
% .verb: display the processing verbally.
% OUTPUT:
% Y: d x n dimension-reduced data matrix.
% ev: 1 - eigenvalues.
% Example: Reduce Swiss roll to 2-D.
% N=1000;
% tt = (3*pi/2)*(1+2*rand(1,N));
% height = 21*rand(1,N);
% X = [tt.*cos(tt); height; tt.*sin(tt)];
% d=2;
% options.k=10;
% [Y,ev] = DRLLE(X, d, options);
%
% CALL: data neighborhood
%
% Algorithm LLE refers to the paper
% S. T. Roweis and L. K. Saul,
% Nonlinear dimensionality reduction
% by locally linear embedding,
        % Initialize inputs.
    options.null=0;
    [~,n] = size(X);
    if ~isfield(options, 'epsilon')
        if isfield(options, 'k')
            K = options.k;
            if numel(K)==1
                K = K*ones(n,1);
            end
        else
            K = 10*ones(n,1);
        end
    end
    if isfield(options, 'nonegativeW')
        nonnegativeW= options.nonnegativeW;
    else
        nonnegativeW=false;
    end
    options.symmetric = 0;
    % Step 1: Construct data-graph.
    fprintf(1,'--Construindo os Grafos\n');
    % find # of neighbors for epsilon-neighborhood.
    [~, nghb] = const_vizinhanca(X,options);
    % Step2: Construct weight matrix
    fprintf(1,'--Construindo a Matrizes de Pesos\n');
    fprintf(1,'If K>D, regularization is used\n');
    W = sparse(n,n);
    if isfield(options, 'epsilon')
    K = sum(nghb, 2);
    end
    tol = 1e-3.*(K > d);
    for ii=1:n
    % shift ith point to origin
    z = X(:,nghb(:,ii))-repmat(X(:,ii),1,K(ii));
    C = z'*z;
    % regularlization (K>D)
    C = C + tol(ii)*trace(C)*eye(K(ii),K(ii));
    W(ii,nghb(:,ii)) = C\ones(K(ii),1);
    end
    if nonnegativeW
    W = W*(W>0);
    end
    % Normorlize W
    W = W./repmat(sum(W,2), 1, n);
    % Step 3: Construct LLE kernel
    fprintf(1,'--Construct LLE kernel.\n');
    K = (speye(n)-W)'*(speye(n)-W);
    % Step 4. Eigen decomposition
    fprintf(1,'--Find dimension-reduced data.\n');
    options.disp = 0; options.isreal = 1;
    options.issym = 1;
    [Y,eigv] = eigs(K,d+1,'SM',options);
    % bottom evect is [1,1,1,1...] with eval 0
    Y = Y(:,1:d)'*sqrt(n);
    D = diag(eigv);
    ev = 1- sqrt(D(1:d));
end
function Y = dr_isomap_isrlivro(X,ndims,options)
    % Isomap DR
    % SYNTAX:
    % Y = Isomap(X,ndims,options);
    % INPUTS:
    % X: a D x N matrix
    % D = dimension, N = #points
    % ndims: Dimension of output data set.
    % OPTIONS:
    % options.epsilon : construct data graph
    % using epsilon-neighbors.
    % options.k : using k-neighbors.
    % options.verb: display the processing verbally.
    %
    % CALL: data neighborhood
    % dijkstra fast.dll or dijkstra
    % Reference:
    % Tenenbaum, de Silva, and Langford,
    % Isomap code -- (c)1998-2000, Josh Tenenbaum
    %
    % Initialize options.
    options.symmetric = 1;

    if isfield(options, 'verb')
        verb =options.verb;
    else
        verb =1;
    end
    % Part 1: Construct data-graph.
    if verb
        disp('- Construindo o Grafo');
    end

    D2 = const_vizinhanca(X,options);
    % Part 2: Compute graph distance.
    if verb
        disp('- Calculando a Matriz de Grafo');
    end
    N = size(D2,1);
    D1 = sqrt(D2);
    DG = symdijkstra(D1, 1:N);
    % Part 3: Generate Isomap kernel
    if verb
    disp('- Gerando O Kernel do Isomap');
    end
    DG = DG.^2;
    GC = -.5*(DG - sum(DG)'*ones(1,N)/N ...
    - ones(N,1)*sum(DG)/N + sum(DG(:))/(N^2));
    % Part 4: Kernel decomposition.
    if verb
        disp('- Projeção em Baixa Dimensão');
    end
    opt.disp = 0; opt.isreal = 1; opt.issym = 1;
    [Y, val] = eigs(GC, ndims, 'LR', opt);
    for i=1:ndims
        Y(:,i) = Y(:,i)*sqrt(val(i,i));
    end
end
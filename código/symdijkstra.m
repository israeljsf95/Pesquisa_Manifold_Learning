function D = symdijkstra(A,pfrm,alld)
% calculate undirected-graph distance.
% It utilities the symmetry for fast computation.
% SYNTAX: D = symdijkstra(A,s,t)
% INPUTS:
% A: n x n symmetric weighted adjacency matrix
% PFROM: node indices
% default []: from all nodes
% ALLD: 1, to all nodes,
% : [] OR 0 (default), to the same set
% OUTPUT:
% D: D(i,j)is distance from i to j
% if ALLD=0, |s| x |s| symmetric matrix
% if ALLD=1, |s| x n matrix, with the most
% left |s| x |s| symmetric sub-matrix, which
% represents distances between the selected
% point set. The order of points in matrix
% is as follows. The selected |s| points
% come first, followed by other points,
% keeping their orders in the original set.
% Copy right note:
% Modified from Michael G. Kays code: Matlog
% Version 1.3 Aug-2000, Copyright (c) 1998-2000.
%%

% Input Error Checking
error(nargchk(1,3,nargin));
[n,cA] = size(A);
if n ~= cA
    error('A must be a square matrix');
    elseif any(any(A < 0))
        error('A must be non-negative');
end

if nargin < 3
    alld = 0;
    if nargin < 2
        pfrm = (1:n);
    end
end

if size(pfrm,1)~=1
    pfrm = (pfrm(:))';
end
% End (Input Error Checking)
% Make the permutation if the selected points
% are not at 1:m.
m=length(pfrm);
if m<n && any(pfrm~=(1:m))
    p=1:n;
    p(pfrm)=[];
    b=[pfrm,p];
    A = A(b,b);
end
% Initialize the output matrix.
if alld
    cD=n;
else
    cD=m;
end

D = zeros(m,cD);
for i = 1:m
%Computing starts all nodes are unlabeled.
    Di = Inf*ones(n+1-i,1); Di(1) = 0;
    NdLbl = 1:(n+1-i);
    NotOut = true(n+1-i,1);
    for k=1:(n+1-i) % Forward phase:
        % compute all path distances
        % through new points
        % Node selection
        [Dj,jj] = min(Di(NotOut));
        j = NdLbl(jj);
        NdLbl(jj) = [];
        NotOut(j) = false;
        [jA,kA,Aj] = find(A(:,j));
        if ~isempty(Aj)
            Dk = Dj + Aj;
            Di(jA) = min(Di(jA),Dk);
        end
    end
    D(i,i:cD) = (Di(1:cD+1-i))';
    %Backward phase: compute path distances
    % through the Out-points
    if i>1
        D(i,i:cD)= min(D(1:i,i:cD)+ ...
        repmat(D(1:i,i),1,cD-i+1));
    end
    %update weighted matrix for
    %avoiding repeated computing.
    A(1,:)=[]; A(:,1)=[];
end
% make the first block symmetric.
D(1:m,1:m)=D(1:m,1:m)+(D(1:m,1:m))';

end
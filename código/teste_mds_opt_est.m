clear all, close all

%Tentativa de Deixar o Ajuste estocastico


%Matriz de distancias D
D = [0	205	3366	2580	3030;...
	205	0	3521	2783	3097;...
 	3366	3521	0	2821	1630;...
	2580	2783	2821	0	3794;...
	3030	3097	1630	3794	0]
nomes = {'aracaju','maceio','rio branco','porto alegre','boa vista'};
% D = d;

% Gerando os pontos P
eta = 10;
d = 2;
p = rand(d, size(D,1));
cont = 2
while (1)
    M = euclid_dist_linalg(p', 2);%matriz de distancias dos pontos
    %embaralhando
    idx = randperm(size(D,1)');
    for k =1:size(D,1)
        nab_J = zeros(d,1);
        v = M(idx(k),:)';
        delta = D(idx(k),:)';
        w = v - delta;
        norm_w = norm(w);
        for i = 1:d
            for j = 1:size(D,2)
               if (j==idx(k))
                   continue
               end
               nab_J(i) = nab_J(i)+ (p(i,idx(k))-p(i,j))*(1-(delta(j)/v(j))); 
            end 
        end
        nab_J = nab_J*1/(norm_w);
        
        p(:,idx(k)) = p(:,idx(k)) - eta*nab_J;
        subplot(211)
        drawnow
        plot(p(1,:),p(2,:),'.')
        text(p(1,:),p(2,:), nomes)
        axis square
    end
%     erro(cont) = sum((D(:)-M(:)).^2)/size(D,1);
    J(cont) = norm_w;
    disp(J(end))
%     if abs(J(end) - J(end-1)) < 1e-3
%         break
%     end
   cont = cont + 1;
   subplot(212)
   drawnow
   plot(J)
end
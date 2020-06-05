
%%
%preludio do mds usando algebra linear
%Matriz de Distancias
clear
D3 = [0,205,3366,2580,3030;...
     205,0,3521,2783,3097;...
     3366,3521,0,2821,1630;...
     2580,2783,2821,0,3794;...
     3030,3097,1630,3794,0];
%supond D seja a matriz do quadrado das distancias
n = size(D, 1);
A = (-1/2)*D.^2;
H = eye(n) - (1/n)*(ones(n,n));
%Matriz de centralização
B = (H*A*H);
[auto_vet, auto_val] = eig(B);
s = struct();
for i=1:n
    a_vet = auto_vet(1:end,1:i);
    a_val = auto_val(1:i,1:i);
    s(i).Y = a_vet*(a_val.^0.5)';
    euc_Y = euclid_dist_linalg(real(s(i).Y),2);
    J(i) = sum(sum((D-euc_Y).^2))/n^2;
end
%J = J/max(J);
subplot(211)
plot([1:n],J)

%ordenando os auto valores pelo tamanho
%auto_valores = sort(auto_valores, 'descend')
%%
capitais_brasil;
D = d;
n = size(D, 1);
A = (-1/2)*D.^2;
H = eye(n) - (1/n)*(ones(n,n));
B = (H*A*H);
[auto_vet, auto_val] = eig(B);
s = struct();
for i=1:n
    a_vet = auto_vet(1:end,1:i);
    a_val = auto_val(1:i,1:i);
    s(i).Y = a_vet*(a_val.^0.5)';
    euc_Y = euclid_dist_linalg(real(s(i).Y),2);
    J(i) = sum(sum((D-euc_Y).^2))/n^2;
end
%J = J/max(J);
subplot(212)
plot([1:n],J);
J = J/max(J);
%%
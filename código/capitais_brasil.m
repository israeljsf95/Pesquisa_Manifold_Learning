clear
cidade(1).nome = 'Aracaju';
cidade(1).coordenadas = [-10.9472 -37.0731];
cidade(2).nome = 'Belém';
cidade(2).coordenadas = [-1.4558 -48.4902];
cidade(3).nome = 'Belo Horizonte';
cidade(3).coordenadas = [-19.9245 -43.9352];
cidade(4).nome = 'Boa Vista';
cidade(4).coordenadas = [2.8235 -60.6758];
cidade(5).nome = 'Brasília';
cidade(5).coordenadas = [-15.7942 -47.8822];
cidade(6).nome = 'Campo Grande';
cidade(6).coordenadas = [-20.4697 -54.6201];
cidade(7).nome = 'Cuiabá';
cidade(7).coordenadas = [-15.6014 -56.0979];
cidade(8).nome = 'Curitiba';
cidade(8).coordenadas = [-25.4244 -49.2654];
cidade(9).nome = 'Florianópolis';
cidade(9).coordenadas = [-27.5949 -48.5482];
cidade(10).nome = 'Fortaleza';
cidade(10).coordenadas = [-3.7319 -38.5267];
cidade(11).nome = 'Goiânia';
cidade(11).coordenadas = [-16.6869 -49.2648];
cidade(12).nome = 'João Pessoa';
cidade(12).coordenadas = [-7.1195 -34.8450];
cidade(13).nome = 'Macapá';
cidade(13).coordenadas = [0.0356 -51.0705];
cidade(14).nome = 'Maceió';
cidade(14).coordenadas = [-9.6498 -35.7089];
cidade(15).nome = 'Manaus';
cidade(15).coordenadas = [-3.1190 -60.0217];
cidade(16).nome = 'Natal';
cidade(16).coordenadas = [-5.7793 -35.2009];
cidade(17).nome = 'Palmas';
cidade(17).coordenadas = [-10.2491 -48.3243];
cidade(18).nome = 'Porto Alegre';
cidade(18).coordenadas = [-30.0346 -51.2177];
cidade(19).nome = 'Porto Velho';
cidade(19).coordenadas = [-8.7612 -63.9004];
cidade(20).nome = 'Recife';
cidade(20).coordenadas = [-8.0476 -34.8770];
cidade(21).nome = 'Rio Branco';
cidade(21).coordenadas = [-9.9754 -67.8249];
cidade(22).nome = 'Rio de Janeiro';
cidade(22).coordenadas = [-22.9068 -43.1729];
cidade(23).nome = 'Salvador';
cidade(23).coordenadas = [-12.9722 -38.5014];
cidade(24).nome = 'São Luís';
cidade(24).coordenadas = [-2.5391 -44.2829];
cidade(25).nome = 'São Paulo';
cidade(25).coordenadas = [-23.5505 -46.6333];
cidade(26).nome = 'Teresina';
cidade(26).coordenadas = [-5.0920 -42.8038];
cidade(27).nome = 'Vitória';
cidade(27).coordenadas = [-20.2976 -40.2958];

po = vertcat(cidade.coordenadas);
d = dist(po');
r = {cidade.nome};


 D = d;

% Gerando os pontos P
eta = 0.1;
d = 2;
p = rand(d, size(D,1));
cont = 2
while (1)
    M = euclid_dist_linalg(p', 2);%matriz de distancias dos pontos
    for k =1:size(D,1)
        nab_J = zeros(d,1);
        v = M(k,:)';
        delta = D(k,:)';
        w = v - delta;
        norm_w = norm(w);
        for i = 1:d
            for j = 1:size(D,2)
               if (j==k)
                   continue
               end
               nab_J(i) = nab_J(i)+ (p(i,k)-p(i,j))*(1-(delta(j)/v(j))); 
            end 
        end
        nab_J = nab_J*1/(norm_w);
        
        p(:,k) = p(:,k) - eta*nab_J;
        
        subplot(211)
        drawnow
        plot(p(1,:),p(2,:),'.')
        text(p(1,:),p(2,:), r)
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
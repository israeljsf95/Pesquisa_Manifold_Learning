

function [H, prob] = calc_prob_linha_i(d_x_i, sigma, i)
    a = sigma;
    prob = exp(-d_x_i*a);
    prob(i) = 0;
%     for j = 1:length(prob)
%         if prob(j) == 0
%             continue
%         end
%         
%         H = H + log(prob(j));
%     end
    H = log(sum(prob)) + a*sum(d_x_i.*prob)./sum(prob);
    prob = prob./sum(prob);
    
end
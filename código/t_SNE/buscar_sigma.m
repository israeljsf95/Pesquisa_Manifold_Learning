
function [prob_final] = buscar_sigma(d_x_i, i, PERP_nom, tol)

    sigma_min = 0;
    sigma_max = inf;
    sigma = 1;
    log_perp = log(PERP_nom);
    [H, prob] = calc_prob_linha_i(d_x_i, sigma, i);
    H_diff = H - log_perp;
    iter = 0;
    limite_superior = 0;
    
    while and(abs(H_diff) > tol, iter < 100)
        
        if H_diff > 0
            if limite_superior == 1
                sigma_min = sigma;
                sigma = (sigma_min + sigma_max)/2;
            else
                sigma_min  = sigma;
                sigma_max = sigma*4;
            end
        else
            sigma_max = sigma;
            sigma = (sigma_min + sigma_max)/2;
            limite_superior = 1;
        end
        
        [H, prob] = calc_prob_linha_i(d_x_i, sigma, i);
        H_diff = H - log_perp;
        iter = iter + 1;
        
    end
    
    prob_final = prob;

end
function [y] = sig_prime(x)
    y_aux = sigmoid(x);
    y = y_aux.*(ones(size(y_aux)) - y_aux);
end
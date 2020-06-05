function y = tanh_prime(x)
    y = ones(size(x))  - tanh(x).^2;
end
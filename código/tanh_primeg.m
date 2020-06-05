function y = tanh_primeg(x)
    y = gpuArray(1  - arrayfun(@tanh, x).^2)w;
end
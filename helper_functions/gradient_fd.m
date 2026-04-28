
function g = gradient_fd(fun, x)

    x = x(:);
    n = numel(x);
    
    f0 = fun(x);
    g = zeros(n,1);
    
    h = 1e-6 * (1 + abs(x));
    
    for i = 1:n
        ei = zeros(n,1);
        ei(i) = 1;
    
        fp = fun(x + h(i)*ei);
        fm = fun(x - h(i)*ei);
    
        g(i) = (fp - fm) / (2*h(i));
    end

end

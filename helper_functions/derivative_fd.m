function df = derivative_fd(fun, s)

    h = 1e-6 * (1 + abs(s));
    
    fp = fun(s + h);
    fm = fun(s - h);
    
    df = (fp - fm) / (2*h);

end


function J = jacobian_fd(fun, x)

    x = x(:);
    n = numel(x);

    f0 = fun(x);
    m = numel(f0);

    J = zeros(m, n);

    h = 1e-6 * (1 + abs(x));

    for i = 1:n
        ei = zeros(n,1);
        ei(i) = 1;

        fp = fun(x + h(i)*ei);
        fm = fun(x - h(i)*ei);

        J(:,i) = (fp - fm) / (2*h(i));
    end

end

% function J = jacobian_fd(f, x)
% 
% fx = f(x);
% 
% n = numel(x);
% eps = 1e-6;
% 
% X = repmat(x,1,n) + eps*eye(n);
% F = zeros(numel(fx), n);
% 
% for i=1:n
%     F(:,i) = f(X(:,i));
% end
% 
% J = (F - fx) / eps;
% 
% end
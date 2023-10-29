function [y, Loss, P] = main(F, NITER, Y, r, X)
v = length(F);
[n, c] = size(F{1});
for i= 1: v
    W{i} = eye(c);
end
a = ones(1, v) / v;
YY = (Y' * Y)^-.5;
d = size(X, 2);
P = diag(ones(1, d) / sqrt(d));

for i = 1 :NITER
    %     Update Y
    lambda1 = -1;
    lambda2 = 0;
    while lambda2 - lambda1 > 1e-5
        for k = 1 : v
            obj(k) = a(k)^2 * norm(Y * YY - F{k} * W{k}, 'fro') ^ 2;
        end
        lambda1 = lambda2;
        XY = Y' * X;
        lambda2 = trace(YY * XY * P * XY' * YY) / sum(obj);
        N = zeros(n, c);
        for k = 1 : v
            N = N + 2 * a(k) ^2 * F{k} * W{k};
        end
        Y = coordinate(Y, lambda2, X * (P^.5), N);
        YY = (Y' * Y)^-.5;
    end

    %     Update a
    for k = 1 : v
        beta(k) = norm(Y * YY - F{k} * W{k}, 'fro') ^ 2;
    end
    a = (2 * beta) .^-1 / sum((2 * beta) .^-1);

    %     Update W

    for k = 1 : v
        [U, ~, V] = svd(F{k}' * Y * YY, 'econ');
        W{k} = U * V';
    end


    %     Update P

    XY = X' * Y;
    p = diag(XY * (Y' * Y)^-1 * XY');
    [~, p_des] = sort(p, 'descend');
    p(p_des(r+1 :end)) = 0;
    P = diag(p ./ norm(p));

    Loss(i) = get_obj(F, P, Y, W, X, a);
    if i > 2 && (Loss(i) - Loss(i - 1)) / Loss(i - 1) < 1e-5
        break
    end
end
[~, y] = max(Y, [], 2);
end
function obj = get_obj(F, P, Y, W, X, a)
XY = X' * Y;
obj1 = trace(P * XY * (Y' * Y)^-1 * XY');
for i = 1 : length(F)
    beta(i) = a(i)^2 * norm(Y * (Y' * Y)^-.5 - F{i} * W{i}, 'fro') ^ 2;
end
obj = obj1 / sum(beta);
end
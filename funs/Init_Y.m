%% max \sum a_i^2 || Y' - F_i * W_i || ^ 2 s.t. a' * 1 = 1
function Y = Init_Y(F)
loss0 = inf;
v = length(F);
[n, c] = size(F{1});
for i= 1: v
    W{i} = eye(c);
end
a = ones(1, v) / v;
for it = 1 :20
    Y = zeros(n, c);
    %     Update Y
    P = a(1)^2 * F{1} * W{1};
    for j = 2 : v
        P = P + a(j)^2 * F{j} * W{j};
    end

    [~, idx] = max(P, [], 2);

    Y = full(sparse(1 : n, idx, 1, n, c));
%     loss1 = get_obj(F, Y, W, a);
%     if loss1 - loss0 > 0
%         sprintf("After Y: %f", loss1 - loss0)
%     end
%     loss0 = loss1;
    %     Update a
    for i = 1 : v
        beta(i) =1 / ( 2 * norm(Y - F{i} * W{i}, 'fro') ^ 2);
    end
    a = beta ./ sum(beta);

%     loss1 = get_obj(F, Y, W, a);
%     if loss1 - loss0 > 0
%         sprintf("After a: %f", loss1 - loss0)
%     end
%     loss0 = loss1;
    %     Update W
    for i = 1 : v
        [U, ~, V] = svd(F{i}' * Y);
        W{i} = U * V';
    end
%     loss1 = get_obj(F, Y, W, a);
%     if loss1 - loss0 > 0
%         sprintf("After W: %f", loss1 - loss0)
%     end
%     loss0 = loss1;
    Loss(it) = get_obj(F, Y, W, a);
    if it > 2 && (Loss(it) - Loss(it - 1)) / Loss(it - 1) < 1e-5
        break
    end
end
end
function obj = get_obj(F, Y, W, a)
for i = 1 : length(F)
    beta(i) = a(i) ^ 2 * norm(Y - F{i} * W{i}, 'fro') ^ 2;
end
obj = sum(beta);
end
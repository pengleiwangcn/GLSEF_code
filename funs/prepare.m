function [F, Y, XX] = prepare(X, k, c)
Graph = generate_graph(X, k);
m = length(Graph);
F = cell(1, m);
A = zeros(size(Graph{1}, 1));
XX = [];
for i=1:m
    graph = Graph{i};
    A = A + graph;
    L = diag(sum(graph)) - graph;
    tmp = eig1(L, c+1, 0, 1);
    tmp(:, 1) = [];
%     tmp = eig1(L, c, 0, 1);
    XX = [XX tmp];
    tmp = diag(sum(tmp.^2, 2) .^-.5) * tmp;
    tmp(isnan(tmp))=0;
    F{i} = tmp;
end
Y = full(finchpp(A, c));
XX = diag(sum(XX.^2, 2).^-.5) * XX;
% XX = (XX - min(XX)) ./ (max(XX) - min (XX));
XX(isnan(XX)) = 0;
end
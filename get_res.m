clc
clear
addpath(genpath('./'))
% load('MSRC-v6.mat')
load('3sources.mat')

v = length(X);
c = length(unique(Y));

% standardize the data
for j=1: v
    colnum=size(X{j},2);
    mole = std(X{j},0,2);
    mole(mole==0) = 1;
    X{j}=(X{j}-mean(X{j},2))./mole;
end

k = 5;
tic

% generate graphs
G = generate_graph(X, k);
[F, Y0, X] = prepare(G, c);

% r is the hyperparameter
for r = c: c: c * v
    tic
    [y, Loss] = main(F, 20, Y0, r, X);
    times = toc;
    res_all(r, :) = [ClusteringMeasure_new(y, Y) times];
end

% display the best result
[~, idx] = max(res_all(:, 2));
res_best = res_all(idx, :)
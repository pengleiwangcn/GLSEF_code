clc
clear

load('MSRC-v6.mat')
% load('3sources.mat')
v = length(X);
for j=1: v
    colnum=size(X{j},2);
    mole = repmat(std(X{j},0,2),1,colnum);
    mole(mole==0) = 1;
    X{j}=(X{j}-repmat(mean(X{j},2),1,colnum))./mole;
end

k = 5;
tic
c = length(unique(Y));
[F, Y0, X] = prepare(X, k, c);

for r = c: c: c * v
    tic
    [y, Loss] = main(F, 20, Y0, r, X);
    times = toc;
    res_all(r, :) = [ClusteringMeasure_new(y, Y) times];
end
[~, idx] = max(res_all(:, 2));
res_best = res_all(idx, :)
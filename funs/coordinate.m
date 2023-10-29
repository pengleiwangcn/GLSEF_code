function Y = coordinate(Y, lambda, B, N)
BtY = B' * Y;
YN = sum(Y .* N);
y = sum(Y);
yBBy = sum(BtY .^ 2, 1);
bb = sum(B .^ 2, 2);
[~, S] = max(Y, [], 2);
S_last = 0;
while any(S_last ~= S)
    S_last = S;
    for i = 1 : size(Y, 1)
        s = S(i);
        if y(s) == 1
            continue
        end
        g1 = (yBBy + 2 * B(i, :) * BtY + bb(i)) ./ (y + 1);
        g2 = yBBy ./ y;
        g1(s) = yBBy(s) / y(s);
        g2(s) = (yBBy(s) - 2* B(i, :) * BtY(:, s) + bb(i)) / (y(s) - 1);

        g3 = (YN + N(i, :)) ./ sqrt(y + 1);
        g4 = YN ./ sqrt(y);
        g3(s) = YN(s) ./ sqrt(y(s));
        g4(s) = (YN(s) - N(i,s)) ./ sqrt(y(s) - 1);
        L = g1 - g2 + lambda * (g3 - g4);
        [~, idx] = max(L);
        if idx ~= s
            Y(i, [s, idx]) = [0, 1];
            S(i) = idx;
            y(s) = y(s) - 1;
            y(idx) = y(idx) + 1;
            BtY(:, s) = BtY(:, s) - B(i, :)';
            BtY(:, idx) = BtY(:, idx) + B(i, :)';
            yBBy(s) = BtY(:, s)' * BtY(:, s);
            yBBy(idx) = BtY(:, idx)' * BtY(:, idx);
            YN(s) = YN(s) - N(i, s);
            YN(idx) = YN(idx) + N(i, idx);
        end
    end
end

end
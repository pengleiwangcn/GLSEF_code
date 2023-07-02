function y_init = finchpp(A0, c)
    warning off
    % graph FINCH initialization
    [esti_y, esti_num_clust] = FINCH(A0, [], 0);
    idx = find(esti_num_clust == c, 1);
    if ~isempty(idx)
        y_init = esti_y(:, idx);
    elseif any(esti_num_clust > c)
        refine_starter = find(esti_num_clust > c, 1, 'last');
        y_init = req_numclust(esti_y(:, refine_starter), A0, c);
    else
        sprintf('FINCH failed to find cluster');
        y_init = kmeans(A0,c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off'); 
    end
    y_init = ind2vec(y_init')';
end

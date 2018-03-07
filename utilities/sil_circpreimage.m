function sel_cands     = sil_circpreimage(silhouette, silh_norm, ...
                                          sigma_sil, gamma, sigma_norm, ...
                                          cand_xys, cand_dists, cand_norms)

% SIL_CIRCPREIMAGE  Finds a circular contour preimage using Appleton and
%   Sun's branch-and-bound method.
                                      
    cs                 = size(cand_xys, 1);
    S                  = size(silhouette, 1);
    start_cands        = [ 1 cs ];
    lower_bounds       = 0;
    circ_val           = Inf;
    
    while ~isempty(lower_bounds) && lower_bounds(1) < circ_val
        [path val] = sil_conspreimage(silhouette, silh_norm, ...
                                      sigma_sil, gamma, sigma_norm, ...
                                      cand_xys, cand_dists, cand_norms, ...
                                      start_cands(1, :), ...
                                      start_cands(1, :), ...
                                      [1 S], true);

        if path(1) == path(S + 1)
            if val < circ_val
                sel_cands = path(1:S);
                circ_val  = val;
            end
        else
            lower_bounds = [ lower_bounds val val ];
            split        = floor((path(1) + path(S + 1)) / 2);
            start_cands  = [ start_cands
                             start_cands(1, 1) split
                             split + 1         start_cands(1, 2) ];
        end
        
        start_cands(1, :) = [];
        lower_bounds(1)   = [];
        
        if path(1) ~= path(S + 1)
            [lower_bounds order] = sort(lower_bounds);
            start_cands          = start_cands(order, :);
        end
    end
    
    %#ok<*AGROW>
end

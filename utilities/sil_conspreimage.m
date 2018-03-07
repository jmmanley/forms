function [path en_val] = sil_conspreimage(silhouette, silh_norm, ...
                                          sigma_sil, gamma, sigma_norm, ...
                                          cand_xys, cand_dists, cand_norms, ...
                                          const_start, const_end, ...
                                          path_fromto, circular)
                                      
% SIL_CONSPREIMAGE  Finds a constrained contour preimage using Dynamic
%   Programming.

cs                 = size(cand_xys, 1);
S                  = size(silhouette, 1);
if path_fromto(2) <= path_fromto(1)
    path_len       = S + path_fromto(2) - path_fromto(1) + 1;
else
    path_len       = path_fromto(2) - path_fromto(1) + 1;
end
path               = zeros(1, path_len + 1);
energy             = zeros(path_len + 1, cs);
previous           = zeros(path_len    , cs);
scaled_cand_dists  = 2 * gamma * (cand_dists .^ 2);
norm_pen1          = cand_norms(:, 3)' / sigma_norm;
norm_pen1          = norm_pen1 .^ 2;

energy(1, :)       = Inf;
energy(1, const_start(1):const_start(2)) = 0;

sz                 = 0;
s                  = path_fromto(1);
done               = false;
while ~done
    sz = sz + 1;
    
    sil_dist   = (silhouette(s .* ones(1, cs), :) - cand_xys) / sigma_sil;
    sil_dist   = sum(sil_dist .^ 2, 2)';
    if s == path_fromto(2)
        sil_dist(1:const_end(1) - 1)  = Inf;
        sil_dist(const_end(2) + 1:cs) = Inf;
    end
    norm_pen2  = (silh_norm(s .* ones(1, cs), 1) - cand_norms(:, 1)) / ...
                 sigma_norm;
    norm_pen2  = norm_pen2' .^ 2;
    norm_pen3  = (silh_norm(s .* ones(1, cs), 2) - cand_norms(:, 2)) / ...
                 sigma_norm;
    norm_pen3  = norm_pen3' .^ 2;

    if circular || sz > 1
        routesofar = energy(sz, :)';
        
        for t = 1:cs
            x = scaled_cand_dists(:, t) + routesofar;
            [energy(sz + 1, t) previous(sz, t)] = min(x);
        end
    end
    energy(sz + 1, :) = energy(sz + 1, :) + ...
                        sil_dist + norm_pen1 + norm_pen2 + norm_pen3;
                    
    if ~circular && sz == 1
        energy(2, :) = energy(2, :) + energy(1, :);
    end
    
    if sz > 1 && s == path_fromto(2)
        done = true;
    else
        s = mod(s, S) + 1;
    end
end

% Build path by looking backwards over minimum-energy paths
[en_val path(path_len + 1)] = min(energy(path_len + 1, :));
for s = path_len:-1:1
    path(s) = previous(s, path(s + 1));
end

if ~circular
    path = path(2:end);
end
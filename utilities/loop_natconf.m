function nc = loop_natconf(valency)

% LOOP_NATCONF  Returns the natural configuration (subdominant eigenvectors)
%   for the Loop subdivision scheme.
%
%   See also LOOP_EVALCHARMAP

% --
% Memoize for performance. All valencies must lie between 3 and max_val
% --
max_val = 20;
persistent nc_store;

if isempty(nc_store)
    nc_store = zeros(max_val + 6, 2, max_val - 2);
end

if nc_store(1, 1, valency - 2) > 0
    nc = nc_store(1:valency + 6, :, valency - 2);
else
    nc        = zeros(valency + 6, 2);
    [~, nzf1] = loop_fouriercomp(valency);
    [v, d]    = eig(nzf1);
    [~, lind] = max(diag(d));
    lev       = v(:, lind) ./ v(1, lind);
    omega     = exp(2 * pi * sqrt(-1) / valency);

    nccomp    = zeros(valency + 6, 1);
    nccomp(1) = lev(2) * omega;
    nccomp(2) = lev(3);
    nccomp(3) = lev(2);
    nccomp(4) = lev(3) * omega;
    nccomp(5) = lev(1) * omega;
    nccomp(6) = lev(1);
    nccomp(7) = lev(3) * omega ^ -1;
    nccomp(8) = lev(1) * omega ^ 2;
    nccomp(9) = 0;

    if valency > 3
        nccomp(10) = lev(1) * omega ^ -1;
    end

    for v = 5:valency;
        nccomp(v + 6) = lev(1) * omega ^ (v - 2);
    end

    nc(:, 1) = real(nccomp);
    nc(:, 2) = imag(nccomp);
    
    nc_store(1:valency + 6, :, valency - 2) = nc;
end

end


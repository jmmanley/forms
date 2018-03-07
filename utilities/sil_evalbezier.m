function pts = sil_evalbezier(pts, t)

% SIL_EVALBEZIER  Evaluate a Bezier curve using the De Casteljau algorithm

degree = size(pts, 1) - 1;

for step = 1:degree
    affine_combs = repmat(1 - t, 1, degree - step + 1);
    smooth_mtx = [diag(affine_combs) zeros(length(affine_combs), 1)];
    smooth_mtx = smooth_mtx + [zeros(length(affine_combs), 1) ...
                               diag(1 - affine_combs)];
    pts = smooth_mtx * pts;
end
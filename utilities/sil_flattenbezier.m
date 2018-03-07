function [ lengths split_pts flat_pts ] = sil_flattenbezier(pts, error)

% SIL_FLATTENBEZIER  Split Bezier curve until it is well-approximated by
%   straight-line segments.

chord = sqrt(sum((pts(4, :) - pts(1, :)) .^ 2));
hull  = sum(sqrt(sum((pts(2:4, :) - pts(1:3, :)) .^ 2, 2)));

if hull - chord < error
    flat_pts  = pts;
    lengths   = (chord + hull) / 2;
    split_pts = [];
else
    left_pts  = zeros(4, 2);
    right_pts = zeros(4, 2);
    
    left_pts(1, :) = pts(1, :);
    left_pts(2, :) = 0.5 * sum(pts(1:2, :));
    left_pts(3, :) = 0.25 * pts(1, :) + 0.5 * pts(2, :) + 0.25 * pts(3, :);
    left_pts(4, :) = 0.125 * pts(1, :) + 0.375 * pts(2, :) + ...
                     0.375 * pts(3, :) + 0.125 * pts(4, :);
                 
    right_pts(1, :) = left_pts(4, :);
    right_pts(2, :) = 0.25 * pts(2, :) + 0.5 * pts(3, :) + 0.25 * pts(4, :);
    right_pts(3, :) = 0.5 * sum(pts(3:4, :));
    right_pts(4, :) = pts(4, :);
    
    [length_left  split_left  flat_left ] = sil_flattenbezier(left_pts,  error);
    [length_right split_right flat_right] = sil_flattenbezier(right_pts, error);
    
    lengths   = [ length_left length_right ];
    split_pts = [ 0.5 * split_left 0.5 0.5 + 0.5 * split_right ];
    flat_pts  = cat(3, flat_left, flat_right);
end

end
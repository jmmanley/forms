function [ x y s ] = align_rotated(targets, rotated_pts)

% ALIGN_ROTATED  Solve least-squares problem on the point constraints.
%
%   See also ALIGN_PROJECT

K = size(targets, 1);
assert(K == size(rotated_pts, 1));
assert(size(targets, 2) == 2);
assert(size(rotated_pts, 2) == 2);

targets = targets';
rotated_pts = rotated_pts';
result = [ repmat(eye(2), K, 1) rotated_pts(:) ] \ targets(:);

x = result(1);
y = result(2);
s = result(3);

end

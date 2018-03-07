function [ center edges ] = loop_limitstencil(valency)

% LOOP_LIMITSTENCIL  Returns stencil (affine weights) for finding a point
%   on a Loop subdivision surface that is the limit of an extraordinary
%   vertex
%
%   See also LIMIT_EXTRAORDINARY

    alpha = (0.375 + cos(2 * pi / valency) / 4) ^ 2 + 0.375;

    subdiv = zeros(valency + 1);
    
    subdiv(1, 1) = alpha;
    subdiv(1, 2:end) = (1 - alpha) / valency;
    subdiv(2:end, 1) = 6;
    subdiv(sub2ind(size(subdiv), 2:valency + 1, 2:valency + 1)) = 6;
    subdiv(sub2ind(size(subdiv), 2:valency + 1, [3:valency + 1 2])) = 2;
    subdiv(sub2ind(size(subdiv), 2:valency + 1, [valency + 1 2:valency])) = 2;
    
    subdiv = subdiv ./ repmat(sum(subdiv, 2), 1, valency + 1);
    
    [v d] = eig(subdiv.');
    [~, limitind] = max(diag(d));
    stencil = v(:, limitind) ./ sum(v(:, limitind));
    
    center = stencil(1);
    edges = sum(stencil(2:end)) / valency;
end


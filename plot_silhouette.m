function plot_silhouette(project, i, res, color)

% PLOT_SILHOUETTE  Plots the continuous silhouette for instance 'i'.

if nargin < 3, res   = 1000; end
if nargin < 4, color = 'r';  end

silhouette = project.images(i).points;
sil_pts    = zeros(res, 2);
sil_params = sil_sample(res, silhouette);

for s = 1:res
    seg              = floor(sil_params(s));
    sil_pts(s, :)    = sil_evalbezier(silhouette(:, :, 1 + seg), ...
                                      sil_params(s) - seg);
end

plot(sil_pts([1:end 1], 1), sil_pts([1:end 1], 2), 'Color', color,'linewidth',1.5);
axis equal;

end
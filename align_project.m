function project = align_project(project,m)

% ALIGN_PROJECT  Ignore the supplied translation and scale for each image,
%   and find these using a least-squares problem on the point constraints
%   instead.
%
% m = a vector containing the indices of images to be aligned

if nargin < 2
    toAlign = 1:length(project.images);
else
    toAlign = m;
end

for ii = 1:length(toAlign)
    i = toAlign(ii);
    image = project.images(i);
    K = length(image.constraints3d);
    
    if K == 0
        continue;
    end
    
    pts = zeros(K, 3);
    for k = 1:K
        pts(k,:) = project.vertices(image.constraints3d(k),:);
    end
    rotated_pts = (image.rotate(1:3, 1:3) * pts')';
    [x y s] = align_rotated(image.constraints2d, rotated_pts(:, 1:2));
    
    if s<0
        disp('scale less than zero!');
        s = abs(s);
    end
    
    image.translate        = eye(4);
    image.translate(1, 4)  = x;
    image.translate(2, 4)  = y;
    image.scale            = s * eye(4);
    image.scale(4, 4)      = 1;
    image.transform        = image.translate * image.rotate * image.scale;
    
    project.images(i) = image;
    
    if isempty(project.images(i).points)
        nPointsPerSilhouette = floor(length(project.images(i).silhouette)/4);
        project.images(i).points = zeros(4,2,nPointsPerSilhouette);
        
        for j=1:4
            idx = [1:4:length(project.images(i).silhouette)-j+1]+j-1;
            idx = idx(1:nPointsPerSilhouette);
            project.images(i).points(j,1,:) = project.images(i).silhouette(idx,1)';
            project.images(i).points(j,2,:) = project.images(i).silhouette(idx,2)';
        end
    end
end

end

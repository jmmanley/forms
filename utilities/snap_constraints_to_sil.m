function project = snap_constraints_to_sil(project,m)

if nargin < 2
    toSnap = [1:length(project.images)];
else
    toSnap = m;
end

for ii=1:length(toSnap)
    i = toSnap(ii);
    
    project.images(i).constraints2d_original = project.images(i).constraints2d;
    
    for j=1:size(project.images(i).constraints2d,1)
        if project.images(i).constraintsonsil(j)
            d = sqrt(sum((project.images(i).constraints2d(j,:) - project.images(i).silhouette).^2,2));
            idx = argmin(d);
            project.images(i).constraints2d(j,:) = project.images(i).silhouette(idx,:);
        end
    end
end

end
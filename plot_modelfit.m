function transformed = plot_modelfit(project, i, model, showsilhouette, plotmesh)

% PLOT_MODELFIT  Plots the fit of the model for instance 'i'.
%   
%   To view the initial fit before optimization, do not pass 'model'.

if nargin < 4
    showsilhouette = true;
end

if nargin < 5
    plotmesh = false;
end

if showsilhouette
    plot_silhouette(project, i); set(gca,'ydir','reverse');
end

P = size(project.vertices, 1);
if nargin < 3 || isempty(model)
    verts      = project.vertices;
    rotscale   = eye(4);
    translate  = eye(4);
else
    rotscale          = eye(4);
    rotscale(1:3,1:3) = model.scale(i) .* model.rotate{i}(1:3,1:3);
    translate         = model.translate{i};
    
    verts      = reshape(model.shapemodes * [1 ; model.shapevars{i}], ...
        P, 3);
end


DT          = translate * ...
    project.images(i).transform * rotscale;
transformed = [ verts ones(P, 1) ] * DT';

if plotmesh
    % Plot mesh outline

    plot_looplimit(project.mesh, transformed(:, 1:3)); set(gca,'ydir','reverse');
    
else
    % Plot 3D surface

    res = 10;
    
    triindex = zeros(res ^ 2, 3);
    currix   = 1;
    
    for x = 0:res - 1
        for y = 0:x
            triindex(currix, :) = [ x * (x + 1) / 2 + y + 1 ...
                (x + 1) * (x + 2) / 2 + y + 1 ...
                (x + 1) * (x + 2) / 2 + y + 2 ];
            currix = currix + 1;
            
            if y < x
                triindex(currix, :) = [ x * (x + 1) / 2 + y + 1 ...
                    (x + 1) * (x + 2) / 2 + y + 2 ...
                    x * (x + 1) / 2 + y + 2 ];
                currix = currix + 1;
            end
        end
    end
    
    tris     = zeros((res + 1) * (res + 2) / 2, 3);
    
    hold on;
    for face = 1:length(project.mesh.edges) / 3
        for x = 0:res
            for y = 0:x
                limitpoint = project.mesh.limitevaluation(face, ...
                    (x - y) / res, y / res);
                currix     = x * (x + 1) / 2 + y + 1;
                tris(currix, :) = limitpoint * transformed(:, 1:3);
            end
        end
        
        trisurf(triindex, tris(:, 1), tris(:, 2), tris(:, 3),'FaceColor',[0.6 0.6 0.7],'EdgeColor','none','AmbientStrength',0.5,'FaceAlpha',0.8);
    end
    hold off; axis equal; set(gca,'ydir','reverse')
    l = light('Position',[0 0 10]);
    l = light('Position',[0 0 -10]);
    lighting gouraud
end

    alpha = 0.08;
    minX = min(transformed(:,1)); maxX = max(transformed(:,1));
    xlim([minX - alpha*abs(minX) maxX + alpha*abs(maxX)]);

    minY = min(transformed(:,2)); maxY = max(transformed(:,2));
    ylim([minY - alpha*abs(minY) maxY + alpha*abs(maxY)]);

    minZ = min(transformed(:,3)); maxZ = max(transformed(:,3));
    zlim([minZ - alpha*abs(minZ) maxZ + alpha*abs(maxZ)]);

end


function transformed = plot_model_fitting(project, i, model, shapevars_history, videoPath)

% Makes a video showing the fitting procedure from
% fit_model_to_image_saveIntermediateSteps.m

if nargin < 5
    showsilhouette = true;
end

if nargin < 6
    plotmesh = false;
end

numFrames = length(shapevars_history);
transformed = zeros(model.P,4,numFrames);

for j=1:numFrames
    rotscale          = eye(4);
    rotscale(1:3,1:3) = model.scale(i) .* model.rotate{i}(1:3,1:3);
    translate         = model.translate{i};
    
    verts      = reshape(model.shapemodes * [1 ; shapevars_history{j}], ...
        model.P, 3);
    
    DT          = translate * ...
        project.images(i).transform * rotscale;
    transformed(:,:,j) = [ verts ones(model.P, 1) ] * DT';
end

transformed = transformed(:,1:3,:);

alpha = 0.08;
minX = min(min(transformed(:,1,:))); maxX = max(max(transformed(:,1,:)));

minY = min(min(transformed(:,2,:))); maxY = max(max(transformed(:,2,:)));

minZ = min(min(transformed(:,3,:))); maxZ = max(max(transformed(:,3,:)));


f = figure('visible','off','position',[1 1 800 300]);
subplot(1,2,1);
imshow(project.images(i).image);
hold on; plot(project.images(i).silhouette(:,1),project.images(i).silhouette(:,2),'r','linewidth',2);
title('Original Frame','fontsize',14); axis off;


    
vw = VideoWriter([videoPath 'fitting_' num2str(i) '.mp4'],'MPEG-4');
vw.FrameRate = 5;
vw.open;

for j=1:numFrames
    a=subplot(1,2,2);
    cla(a);
    title(['Fit #' num2str(j)],'fontsize',14); axis off;
    plot_silhouette(project, i); set(gca,'ydir','reverse'); hold on;
    
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
                tris(currix, :) = limitpoint * transformed(:, 1:3, j);
            end
        end
        
        trisurf(triindex, tris(:, 1), tris(:, 2), tris(:, 3),'FaceColor',[0.6 0.6 0.7],'EdgeColor','none','AmbientStrength',0.5,'FaceAlpha',0.8);
    end
    hold off; axis equal off; set(gca,'ydir','reverse')
    l = light('Position',[0 0 10]);
    l = light('Position',[0 0 -10]);
    lighting gouraud
    
    xlim([minX - alpha*abs(minX) maxX + alpha*abs(maxX)]);
    ylim([minY - alpha*abs(minY) maxY + alpha*abs(maxY)]);
    zlim([minZ - alpha*abs(minZ) maxZ + alpha*abs(maxZ)]);
    
    writeVideo(vw,getframe(f));
end

vw.close;
close(f);


end


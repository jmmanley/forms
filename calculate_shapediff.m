function [shapediff, model_silhouette, minDists, idxMatches] = ...
    calculate_shapediff(project,model,i)

% Provides a potential metric for evaluating the model fit. Sum of
% euclidean distances between original silhouette and silhouette of the 3D
% model.

% Jason Manley, Dec 2017

%%% Sample silhouette.
res = 100;
silhouette = project.images(i).points;
sil_pts    = zeros(res, 2);
sil_params = sil_sample(res, silhouette);

for s = 1:res
    seg              = floor(sil_params(s));
    sil_pts(s, :)    = sil_evalbezier(silhouette(:, :, 1 + seg), ...
                                      sil_params(s) - seg);
end


%%% Find model shape.

rotscale          = eye(4);
rotscale(1:3,1:3) = model.scale(i) .* model.rotate{i}(1:3,1:3);
translate         = model.translate{i};

verts      = reshape(model.shapemodes * [1 ; model.shapevars{i}], ...
    model.P, 3);

DT          = translate * ...
    project.images(i).transform * rotscale;
transformed = [ verts ones(model.P, 1) ] * DT';


%%% Build 3D dolphin model & find silhouette.

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

alltris = [];
alltriindex = [];
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
    
    alltriindex = horzcat(alltriindex,triindex' + length(alltris));
    alltris = horzcat(alltris,tris');
end

% THIS IS SKETCHY! There is a better way to do this if you have the Mapping
% toolbox (which I don't currently have).
% Currently: finding model silhouette by analyzing the image of the 3D
% model, with some careful indexing of axis limits.

f = figure('visible','off');
trisurf(alltriindex', alltris(1,:), alltris(2,:), alltris(3,:),'FaceColor',[0 0 0],'EdgeColor','none','AmbientStrength',0.5,'FaceAlpha',1);
set(gca,'position',[0 0 1 1],'units','normalized')
view(2); axis off; grid off; set(gca,'ydir','reverse')
xl = xlim; yl = ylim;

F = getframe(f);
[X] = frame2im(F);
close(f);

mapx = linspace(xl(1),xl(2),size(X,2));
mapy = linspace(yl(1),yl(2),size(X,1));
[px, py] = find(del2(double(rgb2gray(X)))>0);

res = 300;
idx = round(linspace(1,length(px),res));
model_silhouette = zeros(res,2);
model_silhouette(:,2) = mapy(px(idx));
model_silhouette(:,1) = mapx(py(idx));


%%% Calculate difference between model shape and silhouette shape.

% First, find correspondence between model and actual silhouette by
% minimizing distance between points.

idxMatches = zeros(1,length(sil_pts));
minDists = zeros(1,length(sil_pts));
for ii = 1:length(sil_pts)
    dists = sqrt(sum((sil_pts(ii,:) - model_silhouette).^2,2));
    [minDists(ii), idxMatches(ii)] = min(dists);
end

shapediff = sum(minDists);

end






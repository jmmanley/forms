function [movie] = plot_basisshape_sweep(project,model,mode,saveFile,res,low,high,n)

% Makes a movie "sweeping" through values of one basis shape.
%
% Jason Manley, Nov 2017

if nargin < 5
    res = 30;
end

if nargin < 6
    low  = min(combineCells(model.shapevars));
    high = max(combineCells(model.shapevars));
    n = 10;
end

disp('generating plots...')

xx = linspace(low,high,res);
%xx = [xx fliplr(xx)];

shapevars = ones(model.M+1,1);
shapemodes = model.shapemodes;
shapemode  = reshape(shapemodes(:, mode + 1), model.P, 3);

trires = 10;

triindex = zeros(trires ^ 2, 3);    
currix   = 1;

for x = 0:trires - 1
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

for i=1:length(xx)
    shapevars(mode+1) = xx(i);
    transformed{i} = reshape(shapemodes * shapevars,model.P,3)';
end
allTransformed = combineCells(transformed,2);
[mins] = min(allTransformed,[],2); minx = mins(1); miny = mins(2); minz = mins(3);
[maxes] = max(allTransformed,[],2); maxx = maxes(1); maxy = maxes(2); maxz = maxes(3);

figure('Position',[800 500 1000 400],'Visible','Off');
subplot(1,2,1);
set(gca,'ydir','reverse')
l = light('Position',[0 0 10]);
l = light('Position',[0 0 -10]);
lighting gouraud

subplot(1,2,2);
set(gca,'ydir','reverse')

for i=1:length(xx)
    tris     = zeros((trires + 1) * (trires + 2) / 2, 3);
    modesize = zeros((trires + 1) * (trires + 2) / 2, 1);
    
    subplot(1,2,1); cla; hold on; subplot(1,2,2); cla; hold on;
    for face = 1:length(project.mesh.edges) / 3
        for x = 0:trires
            for y = 0:x
                limitpoint = project.mesh.limitevaluation(face, ...
                    (x - y) / trires, y / trires);
                currix     = x * (x + 1) / 2 + y + 1;
                
                tris(currix, :) = limitpoint * transformed{i}';
                
                modevec = limitpoint * shapemode;
                modesize(currix) = norm(modevec, 2);
            end
        end
        
        subplot(1,2,1);
        trisurf(triindex, tris(:, 1), tris(:, 2), tris(:, 3),'FaceColor',[0.6 0.6 0.7],'EdgeColor','none')
        subplot(1,2,2);
        trisurf(triindex, tris(:, 1), tris(:, 2), tris(:, 3),'FaceColor',[0.6 0.6 0.7],'EdgeColor','none')
    end
    subplot(1,2,1); view(2); axis equal off; %xlim([minx maxx]); ylim([miny maxy]); zlim([minz maxz]);
    l = light('Position',[0 0 10]);
    l = light('Position',[0 0 -10]);
    lighting gouraud
    subplot(1,2,2); view(3); axis equal off; %xlim([minx maxx]); ylim([miny maxy]); zlim([minz maxz]);
    l = light('Position',[0 0 10]);
    l = light('Position',[0 0 -10]);
    lighting gouraud
    movie(i) = getframe(gcf);
end

for i=1:length(xx)-1
    movie(length(xx)+i) = movie(length(xx)-i);
end

movie = repmat(movie,[1 n]);

disp('writing to movie...')

if strcmp(saveFile(end-3:end),'.mp4')
    type = 'MPEG-4';
else
    type = 'Motion JPEG AVI';
end

vw = VideoWriter(saveFile,type);
vw.open;
for i=1:length(movie)
    writeVideo(vw,movie(i));
end
vw.close;

end
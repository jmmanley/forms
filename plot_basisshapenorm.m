function plot_basisshapenorm(project, model, mode)

% PLOT_BASISSHAPENORM  Shows the norm of a basis shape as a coloured
%   version of the zeroth shape mode.

assert(mode > 0);

shapemodes = model.shapemodes;
shapemode  = reshape(shapemodes(:, mode + 1), model.P, 3);
meansurf   = reshape(shapemodes(:, 1), model.P, 3);

res = 8;

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
modesize = zeros((res + 1) * (res + 2) / 2, 1);

hold on;
for face = 1:length(project.mesh.edges) / 3
    for x = 0:res
        for y = 0:x
            limitpoint = project.mesh.limitevaluation(face, ...
                                                      (x - y) / res, y / res);
            currix     = x * (x + 1) / 2 + y + 1;
            tris(currix, :) = limitpoint * meansurf;
            
            modevec = limitpoint * shapemode;
            modesize(currix) = norm(modevec, 2);
        end
    end

    trisurf(triindex, tris(:, 1), tris(:, 2), tris(:, 3), modesize);
end
hold off;
axis equal;
shading interp;

end
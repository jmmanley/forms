function plot_looplimit(mesh, vertices, res, InteriorColor, EdgeColor)

% PLOT_LOOPLIMIT  Draw the limit of a loop surface given by the MESHTRI
%   mesh and the vertices

    if nargin < 3
        res = 2;
    end
    if nargin < 4
        InteriorColor = [ 0.8 0.8 0.8 ];
    end
    if nargin < 5
        EdgeColor     = [ 1   0   0   ];
    end
    
    currcol = zeros(res + 1, 3);
    
    for face = 1:length(mesh.edges) / 3;
        for x = 0:res
            for y = 0:x
                limitpoint = mesh.limitevaluation(face, ...
                                                  (x - y) / res, y / res);
                currcol(y + 1, :) = limitpoint * vertices;
            end
            
            if x > 0
                l = [ prevcol(1, :) ; currcol(1, :) ];
                line(l(:, 1), l(:, 2), l(:, 3), 'Color', EdgeColor);
                
                l = [ prevcol(x, :) ; currcol(x + 1, :) ];
                line(l(:, 1), l(:, 2), l(:, 3), 'Color', EdgeColor);
            end
            
            if ~all(InteriorColor == 1)
                for yv = 1:x - 1
                    l = [ prevcol(yv + 1, :) ; currcol(yv + 1, :) ];
                    line(l(:, 1), l(:, 2), l(:, 3), 'Color', InteriorColor);
                    l = [ prevcol(yv, :) ; currcol(yv + 1, :) ];
                    line(l(:, 1), l(:, 2), l(:, 3), 'Color', InteriorColor);
                end
            end
            
            if x < res
                if ~all(InteriorColor == 1)
                    line(currcol(1:x + 1, 1), ...
                         currcol(1:x + 1, 2), currcol(1:x + 1, 3), ...
                         'Color', InteriorColor);
                end
            else
                line(currcol(:, 1), currcol(:, 2), currcol(:, 3), ...
                     'Color', EdgeColor);
            end
            
            prevcol = currcol;
        end
    end
    axis equal;
end


classdef meshtri < handle
    % MESHTRI  Class used to represent a triangle mesh.
    
    properties
        origverts
        edges
        verts
        wholeedges
    end
    
    methods
        function obj = meshtri(vertices, faces)
            obj.origverts = vertices;
           
            if nargin > 1
                obj.verts = repmat(meshtri_vertex, 1, size(vertices, 1));
                for v = 1:size(vertices, 1)
                    obj.verts(v) = meshtri_vertex(size(vertices, 1), v);
                end

                pairtable      = zeros(3 * size(faces, 1), 2);
                obj.edges      = repmat(meshtri_halfedge, ...
                                        1, 3 * size(faces, 1));
                for f = 1:size(faces, 1)
                    % Every third edge defines a face.
                    edge1 = meshtri_halfedge(obj.verts(faces(f, 1)));
                    edge2 = meshtri_halfedge(obj.verts(faces(f, 2)));
                    edge3 = meshtri_halfedge(obj.verts(faces(f, 3)));

                    edge1.next = edge2;
                    edge2.next = edge3;
                    edge3.next = edge1;

                    obj.edges(3 * f - 2) = edge1;
                    obj.edges(3 * f - 1) = edge2;
                    obj.edges(3 * f    ) = edge3;
                    
                    edge1.index_in_mesh  = 3 * f - 2;
                    edge2.index_in_mesh  = 3 * f - 1;
                    edge3.index_in_mesh  = 3 * f;

                    pairtable(3 * f - 2, :) = [ faces(f, 3) faces(f, 1) ];
                    pairtable(3 * f - 1, :) = [ faces(f, 1) faces(f, 2) ];
                    pairtable(3 * f    , :) = [ faces(f, 2) faces(f, 3) ];
                end
                
                obj.wholeedges = repmat(meshtri_halfedge, ...
                                        1, 3 * size(faces, 1) / 2);
                wholeedgeix    = 1;
                for p = 1:size(pairtable, 1)
                    obj.edges(p).pair = ...
                        obj.edges(pairtable(:, 1) == pairtable(p, 2) & ...
                                  pairtable(:, 2) == pairtable(p, 1));

                    if isempty(obj.edges(p).edge)
                        obj.edges(p).edge = meshtri_edge;
                        obj.edges(p).pair.edge = obj.edges(p).edge;
                        obj.wholeedges(wholeedgeix) = obj.edges(p);
                        wholeedgeix = wholeedgeix + 1;
                    end

                    if isempty(obj.edges(p).vert.edge)
                        obj.edges(p).vert.edge = obj.edges(p).pair;
                    end
                end
                
                for v = 1:size(vertices, 1)
                    obj.verts(v).setValency();
                end
            end
        end
        
        function energy = thinplate(obj)
            energy = zeros(size(obj.origverts, 1));
            
            for e = 3:3:length(obj.edges)
                % Middle triangle
                defs = meshtri.middle_triangle(obj.edges(e));
                energy = energy + defs' * energy_ordinary() * defs;
                
                % Edge triangles
                p = obj.edges(e);
                for i = 1:3
                    [defs valency] = meshtri.edge_triangle(p);
                              
                    if valency == 6
                        energy = energy + ...
                                 defs' * energy_ordinary() * defs;
                    else
                        energy = energy + defs' * ...
                                 energy_extraordinary(valency) * defs;
                    end
                         
                    p = p.next;
                end
            end
    
            % Lost symmetry due to round-off error?
            energy = (energy + energy') / 2;
        end
        
        function [def deriv sec] = limitevaluation(obj, face, u, v)            
            if isa(face, 'meshtri_halfedge')
                topedge = face;
            else
                topedge = obj.edges(face * 3);
            end
            
            if v <= 0.5 && u <= 0.5 && u + v >= 0.5
                defs = meshtri.middle_triangle(topedge);
                if nargout > 2
                    [def deriv sec] = limit_ordinary(1 - 2 * u, ...
                                                     1 - 2 * v);
                else
                    [def deriv]     = limit_ordinary(1 - 2 * u, ...
                                                     1 - 2 * v);
                end
                deriv            = -deriv;
            else
               if u > 0.5
                    vertedge = topedge;
                    u   = 2 * u - 1;
                    v   = v * 2;
                    rot = 2;
                    
                    nu  = v;
                    nv  = 1 - (u + v);
                elseif v > 0.5
                    vertedge = topedge.next.next;
                    u   = u * 2;
                    v   = 2 * v - 1;
                    rot = 1;
                    
                    nu  = 1 - (u + v);
                    nv  = u;
                elseif u + v < 0.5
                    vertedge = topedge.next;
                    nu  = u * 2;
                    nv  = v * 2;
                    rot = 0;
                end
                
                [defs valency] = meshtri.edge_triangle(vertedge);
                
                if nargout > 2
                    if valency == 6
                        [def deriv sec] = limit_ordinary(nu, nv);
                    else
                        [def deriv sec] = limit_extraordinary(nu, nv, ...
                                                              valency);
                    end
                else
                    if valency == 6
                        [def deriv]     = limit_ordinary(nu, nv);
                    else
                        [def deriv]     = limit_extraordinary(nu, nv, ...
                                                              valency);
                    end
                end
                
                if rot == 1
                    deriv = [ deriv(2, :) ; deriv(3, :) ; deriv(1, :) ];
                    if nargout > 2
                        sec   = [ sec(3, :)                             ;
                                  sec(3, :) - sec(2, :)                 ;
                                  sec(3, :) - 2 * sec(2, :) + sec(1, :) ];
                    end
                elseif rot == 2
                    deriv = [ deriv(3, :) ; deriv(1, :) ; deriv(2, :) ];
                    if nargout > 2
                        sec   = [ sec(1, :) - 2 * sec(2, :) + sec(3, :) ;
                                  sec(1, :) - sec(2, :)                 ;
                                  sec(1, :)                             ];
                    end
                end
            end
            
            def     = def   * defs;
            deriv   = deriv * defs * 2;
            if nargout > 2
                sec = sec   * defs * 4;
            end
        end
        
        function allocNewVerts(obj)
            if isempty(obj.verts(1).newVert)
                for v = 1:length(obj.verts)
                    obj.verts(v).setNewVert();
                end
                for e = 1:length(obj.edges)
                    obj.edges(e).setNewVert();
                end
            end
        end
    end
    
    methods(Static)
        function defs = middle_triangle(e)
            defs = [ e.next.pair.next.next.edge.newVert.definition ;
                     e.next.vert.newVert.definition ;
                     e.next.next.pair.next.edge.newVert.definition ;
                     e.next.pair.next.edge.newVert.definition ;
                     e.next.edge.newVert.definition ;
                     e.next.next.edge.newVert.definition ;
                     e.next.next.pair.next.next.edge.newVert.definition ;
                     e.vert.newVert.definition ;
                     e.edge.newVert.definition ;
                     e.pair.vert.newVert.definition ;
                     e.pair.next.next.edge.newVert.definition ;
                     e.pair.next.edge.newVert.definition ];
        end
        
        function [defs valency] = edge_triangle(e)
            valency = e.vert.valency;
            n       = e.next;
            
            defs = [ n.vert.newVert.definition ;
                     n.next.edge.newVert.definition ;
                     n.next.vert.newVert.definition ;
                     n.pair.next.next.edge.newVert.definition ;
                     n.edge.newVert.definition ;
                     e.edge.newVert.definition ;
                     e.pair.next.edge.newVert.definition ;
                     n.pair.next.edge.newVert.definition ;
                     e.vert.newVert.definition ;
                     e.pair.next.next.edge.newVert.definition ];

            if valency == 3
                defs      = defs(1:9, :);
            elseif valency > 4
                defs      = [ defs ; ...
                              zeros(valency - 4, size(defs, 2)) ];
                
                start = n.pair.next.pair;
                for v = 5:valency
                    start          = start.next.pair;
                    defs(v + 6, :) = start.edge.newVert.definition;
                end
            end
        end
        
        function [dist du1 dv1 du2 dv2] = distbetween(edge1, u1, v1, ...
                                                      edge2, u2, v2)
            w1   = 1 - (u1 + v1);
            w2   = 1 - (u2 + v2);
            rot1 = 0;
            rot2 = 0;
            
            % Search current and adjacent triangles of edge1 for edge2
            if edge1.isInLoop(edge2)
                % Same triangle
                if edge1.next == edge2
                    u2   = v2;
                    v2   = w2;
                    rot2 = 2;
                elseif edge1.next.next == edge2
                    v2   = u2;
                    u2   = w2;
                    rot2 = 1;
                end
                
                [dist du1 dv1 du2 dv2] = ...
                    meshtri.dist_same_tri(u1, v1, u2, v2, ...
                                          [edge1.vert.valency ...
                                           edge1.next.vert.valency ...
                                           edge1.next.next.vert.valency]);
            else
                found_edge = false;
                
                if edge1.pair.isInLoop(edge2)
                    edge1 = edge1.next.next;
                    u1    = v1;
                    v1    = w1;
                    rot1  = 2;
                    found_edge = true;
                elseif edge1.next.pair.isInLoop(edge2)
                    found_edge = true;
                elseif edge1.next.next.pair.isInLoop(edge2)
                    edge1 = edge1.next;
                    v1    = u1;
                    u1    = w1;
                    rot1  = 1;
                    found_edge = true;
                end
                
                if found_edge
                    % edge1.next.pair is in the same triangle as edge2
                    
                    if edge1.next.pair == edge2
                        edge2 = edge2.next;
                        v2    = u2;
                        u2    = w2;
                        rot2  = 1;
                    elseif edge1.next.pair == edge2.next
                        edge2 = edge2.next.next;
                        u2    = v2;
                        v2    = w2;
                        rot2  = 2;
                    end
                    
                    assert(edge1.next.pair.next == edge2);
                    t = v1 / (u2 + v1);
                    if isnan(t)
                        t = 0;
                    end
                    
                    edge1_vals = [ edge1.vert.valency
                                   edge1.next.vert.valency
                                   edge1.next.next.vert.valency ];
                    edge2_vals = [ edge2.vert.valency
                                   edge2.next.vert.valency
                                   edge2.next.next.vert.valency ];
                    
                    [dist1 du11 dv11 du21] = ...
                        meshtri.dist_same_tri(u1, v1, (1-t) * u1 + t * v2, ...
                                              0, edge1_vals);
                    [dist2, ~, dv12 du22 dv22] = ...
                        meshtri.dist_same_tri(0, (1 - t) * u1 + t * v2, ...
                                              u2, v2, edge2_vals);
                    dist = dist1 + dist2;
                    du1  = du11 + (du21 + dv12) * (1 - t);
                    dv1  = dv11 + (du21 + dv12) * u2 * (v2 - u1) / ...
                                                  ((u2 + v1) ^ 2);
                    du2  = du22 + (du21 + dv12) * v1 * (u1 - v2) / ...
                                                  ((u2 + v1) ^ 2);
                    dv2  = dv22 + (du21 + dv12) * t;
                else
                    % Same vertex area?
                    
                    [~, i] = max([u1 w1 v1]);
                    switch i
                        case 1
                            edge1 = edge1.next.next;
                            u1    = v1;
                            v1    = w1;
                            rot1  = 2;
                        case 3
                            edge1 = edge1.next;
                            v1    = u1;
                            u1    = w1;
                            rot1  = 1;
                    end
                    [~, j] = max([u2 w2 v2]);
                    switch j
                        case 1
                            edge2 = edge2.next.next;
                            u2    = v2;
                            v2    = w2;
                            rot2  = 2;
                        case 3
                            edge2 = edge2.next;
                            v2    = u2;
                            u2    = w2;
                            rot2  = 1;
                    end
                    
                    test_edge = edge1.next.next.pair.next;
                    rot_steps = 1;
                    while test_edge ~= edge1.next.pair
                        if test_edge.isInLoop(edge2)
                            break;
                        end
                        
                        rot_steps = rot_steps + 1;
                        test_edge = test_edge.pair.next;
                    end
                    
                    if test_edge.next == edge2
                        valency   = edge1.next.vert.valency;
                        assert(rot_steps > 1 && rot_steps < valency - 1);
                        
                        [p1 d1] = loop_evalcharmap(valency, u1, v1);
                        [p2 d2] = loop_evalcharmap(valency, u2, v2);
                        p2      = complex(p2(1), p2(2));
                        p2      = p2 * exp(complex(0, ...
                                           2 * pi * rot_steps / valency));
                        p2      = [ real(p2) imag(p2) ];
                        d2      = complex(d2(:, 1), d2(:, 2));
                        d2      = d2 * exp(complex(0, ...
                                           2 * pi * rot_steps / valency));
                        d2      = [ real(d2) imag(d2) ];
                        dist    = sum((p1 - p2) .^ 2);
                        du1     = sum(2 * (p1 - p2) .*  d1(2, :));
                        dv1     = sum(2 * (p1 - p2) .* -d1(1, :));
                        du2     = sum(2 * (p2 - p1) .*  d2(2, :));
                        dv2     = sum(2 * (p2 - p1) .* -d2(1, :));
                
                        [dist du1 dv1 du2 dv2] = ...
                            meshtri.sqrt_approx(dist, du1, dv1, du2, dv2);
                    else
                        dist = Inf;
                        du1  = 0;
                        dv1  = 0;
                        du2  = 0;
                        dv2  = 0;
                    end
                end
            end
            
            switch rot1
                case 1
                    tmp = dv1;
                    dv1 = -du1;
                    du1 = tmp - du1;
                case 2
                    tmp = du1;
                    du1 = -dv1;
                    dv1 = tmp - dv1;
            end
            
            switch rot2
                case 1
                    tmp = dv2;
                    dv2 = -du2;
                    du2 = tmp - du2;
                case 2
                    tmp = du2;
                    du2 = -dv2;
                    dv2 = tmp - dv2;
            end
            
            if dist > sqrt(3) / 4
                dist = sqrt(3) / 4; du1 = 0; dv1 = 0; du2 = 0; dv2 = 0;
            end
            if isnan(du1), du1 = 0; end
            if isnan(dv1), dv1 = 0; end
            if isnan(du2), du2 = 0; end
            if isnan(dv2), dv2 = 0; end
        end
        
        function [dist du1 dv1 du2 dv2] = dist_same_tri(u1, v1, u2, v2, ...
                                                        valencies, i, j)
            w1 = 1 - (u1 + v1);
            w2 = 1 - (u2 + v2);
            
            if nargin < 7
                % Which vertex area (u1, v1) is in:
                %   i = 1 : closest to u1 vertex
                %   i = 2 : closest to w1 vertex
                %   i = 3 : closest to v1 vertex
                [~, i] = max([u1 w1 v1]);
                % j is the same but for (u2, v2):
                [~, j] = max([u2 w2 v2]);
            end
            
            if i ~= 2
                switch i
                    case 1
                        u1   = v1;
                        v1   = w1;
                        j    = mod(j, 3) + 1;
                        u2   = v2;
                        v2   = w2;
                        oldi = 1;
                    case 3
                        v1   = u1;
                        u1   = w1;
                        j    = mod(j + 1, 3) + 1;
                        v2   = u2;
                        u2   = w2;
                        oldi = 3;
                end
                w1 = 1 - (u1 + v1);
                w2 = 1 - (u2 + v2);
                i  = 2;
            else
                oldi = 2;
            end
            
            if i == j
                n       = valencies(i);
                [p1 d1] = loop_evalcharmap(n, u1, v1);
                [p2 d2] = loop_evalcharmap(n, u2, v2);
                dist    = sum((p1 - p2) .^ 2);
                du1     = sum(2 * (p1 - p2) .*  d1(2, :));
                dv1     = sum(2 * (p1 - p2) .* -d1(1, :));
                du2     = sum(2 * (p2 - p1) .*  d2(2, :));
                dv2     = sum(2 * (p2 - p1) .* -d2(1, :));
                
                [dist du1 dv1 du2 dv2] = ...
                    meshtri.sqrt_approx(dist, du1, dv1, du2, dv2);
            else
                if abs(u1 - 1/3) < 1e-8 && abs(v1 - 1/3) < 1e-8
                    % The method below fails if a point is directly in the 
                    % middle of a triangle, as then the line through the
                    % centroid is ill-defined. So move the point very
                    % slightly off the centroid.
                    u1 = u1 - 1e-6;
                    v1 = v1 - 1e-6;
                end
                ang1 = cart2pol(u1 - 1/3, v1 - 1/3);
                ang2 = cart2pol(u2 - 1/3, v2 - 1/3);
                t_vw = (w1 - v1) / ((v2 - w2) + (w1 - v1));
                t_wu = (u1 - w1) / ((w2 - u2) + (u1 - w1));
                
                if (mod(ang1 - ang2, 2 * pi) < pi && t_vw >= 0 ...
                    && t_vw <= 1) || isnan(t_wu)
                    newi = 3;
                    isct = t_vw;
                else
                    newi = 1;
                    isct = t_wu;
                end
                
                if isct < 0 || isct > 1
                    % Only if we're calculating a distance which is
                    % directly in line with a vertex area boundary.
                    isct = 0;
                end
                
                [dist1 du11 dv11 du21 dv21] = ...
                    meshtri.dist_same_tri(u1, v1,...
                                          isct * u2 + (1 - isct) * u1, ...
                                          isct * v2 + (1 - isct) * v1, ...
                                          valencies, i, i);
                [dist2 du12 dv12 du22 dv22] = ...
                    meshtri.dist_same_tri(isct * u2 + (1 - isct) * u1, ...
                                          isct * v2 + (1 - isct) * v1, ...
                                          u2, v2, valencies, newi, j);

                dist = dist1 + dist2;
                
                if newi == 3
                    denom = (u1 - u2 + 2 * v1 - 2 * v2) ^ 2;
                    dvdu1 = (v1 - v2) * (u2 + 2 * v2 - 1) / denom;
                    dudu1 = -2 * dvdu1;
                    dvdv1 = (u2 - u1) * (u2 + 2 * v2 - 1) / denom;
                    dudv1 = -2 * dvdv1;
                    dvdu2 = (v2 - v1) * (u1 + 2 * v1 - 1) / denom;
                    dudu2 = -2 * dvdu2;
                    dvdv2 = (u1 - u2) * (u1 + 2 * v1 - 1) / denom;
                    dudv2 = -2 * dvdv2;
                else % newi == 1
                    denom = (2 * u1 - 2 * u2 + v1 - v2) ^ 2;
                    dudu1 = (v2 - v1) * (2 * u2 + v2 - 1) / denom;
                    dvdu1 = -2 * dudu1;
                    dudv1 = (u1 - u2) * (2 * u2 + v2 - 1) / denom;
                    dvdv1 = -2 * dudv1;
                    dudu2 = (v1 - v2) * (2 * u1 + v1 - 1) / denom;
                    dvdu2 = -2 * dudu2;
                    dudv2 = (u2 - u1) * (2 * u1 + v1 - 1) / denom;
                    dvdv2 = -2 * dudv2;
                end
                
                du1 = du11 + (du21 + du12) * dudu1 + (dv21 + dv12) * dvdu1;
                dv1 = dv11 + (du21 + du12) * dudv1 + (dv21 + dv12) * dvdv1;
                du2 = du22 + (du21 + du12) * dudu2 + (dv21 + dv12) * dvdu2;
                dv2 = dv22 + (du21 + du12) * dudv2 + (dv21 + dv12) * dvdv2;
            end
            
            switch oldi
                case 1
                    tmp = du1;
                    du1 = -dv1;
                    dv1 = tmp - dv1;
                    tmp = du2;
                    du2 = -dv2;
                    dv2 = tmp - dv2;
                case 3
                    tmp = dv1;
                    dv1 = -du1;
                    du1 = tmp - du1;
                    tmp = dv2;
                    dv2 = -du2;
                    du2 = tmp - du2;
            end
        end
        
        function [dist du1 dv1 du2 dv2] = sqrt_approx(dist, du1, dv1, ...
                                                            du2, dv2)
            % Approximation to the square root function with finite
            % derivatives
            a       = 1e-1;
            if dist > a^2
                dist = sqrt(dist);
                du1  = du1 / (2 * dist);
                dv1  = dv1 / (2 * dist);
                du2  = du2 / (2 * dist);
                dv2  = dv2 / (2 * dist);
            else
                dndd =    1 / a + 2 * dist / (2 * a^3) ...
                    - 3 * dist^2 / (2 * a^5);
                dist = dist / a +   dist^2 / (2 * a^3) ...
                    -     dist^3 / (2 * a^5);
                du1  = du1 * dndd;
                dv1  = dv1 * dndd;
                du2  = du2 * dndd;
                dv2  = dv2 * dndd;
            end
        end
    end
end


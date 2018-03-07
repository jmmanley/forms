function [ cand_ixs cand_uvs cand_limits ...
           cand_dists cand_derivs ] = init_candidates(mesh, res)

% INIT_CANDIDATES Set up discrete candidate positions for global search
%   for contour preimages.
%
%   See also READ_PROJECT
       
    if nargin < 2, res = 4; end;

    P = length(mesh.verts);
    F = length(mesh.edges) / 3;
    
    no_cands    = P + F * res * (res + 1) / 2;
    cand_ixs    = zeros(no_cands, 1);
    cand_limits = zeros(no_cands, P);
    cand_derivs = zeros(no_cands, P, 2);
    cand_dists  =  ones(no_cands);
    faces_done  = zeros(F);
    
    for p = 1:P
        cand_ixs(p)  = mesh.verts(p).edge.next.index_in_mesh;
    end

    face_edge_table = zeros(F, 1);
    ptsperface = res * (res + 1) / 2;
    for f = 1:F
        face_edge_table(f) = mesh.edges(f * 3).index_in_mesh;
        cand_ixs(P + (f - 1) * ptsperface + 1:P + f * ptsperface) = ...
            face_edge_table(f);
    end
    
    face_uvs   = zeros(ptsperface, 2);
    f = 1;
    for u = 1:res
        for v = 1:res + 1 - u
            face_uvs(f, :) = [ (2 * u - 1) (2 * v - 1) ] / (2 * res + 1);
            
            if 3 * (2 * u - 1) == 2 * res + 1 && ...
               3 * (2 * v - 1) == 2 * res + 1
                face_uvs(f, :) = face_uvs(f, :) - 1e-4;
            end
            f = f + 1;
        end
    end
    cand_uvs   = [zeros(P, 2) ; repmat(face_uvs, F, 1)];
    
    fprintf(1, '\nInitialising candidates in %d faces:\n', F);
    for p = 1:P
        [cand_limits(p, :) deriv] = ...
            mesh.limitevaluation(mesh.edges(cand_ixs(p)), 0, 0);
        cand_derivs(p, :, 1) = deriv(2, :);
        cand_derivs(p, :, 2) = [ -1/sqrt(3) 0 1/sqrt(3) ] * deriv;
        
        e = mesh.edges(cand_ixs(p)).next.next;
        n = mesh.edges(cand_ixs(p)).next.vert.valency;
        for v = 1:n
            face = findedge(e);
            e    = e.pair.next;
            
            for c = 1:ptsperface
                cd  = P + (face - 1) * ptsperface + c;
                cand_dists(p, cd) = ...
                    mesh.distbetween(mesh.edges(cand_ixs(p)), 0, 0, ...
                                     mesh.edges(cand_ixs(cd)), ...
                                     cand_uvs(cd, 1), cand_uvs(cd, 2));
                cand_dists(cd, p) = cand_dists(p, cd);
            end
        end
    end
    
    for f = 1:F
        fprintf(1, '%3d  ', f);
        if ~mod(f, 10), fprintf(1, '\n'); end
        
        % Find the one-ring faces of this face
        v1 = mesh.edges(f * 3).vert.valency;
        v2 = mesh.edges(f * 3).next.vert.valency;
        v3 = mesh.edges(f * 3).next.next.vert.valency;
        onering = zeros(1, v1 + v2 + v3 - 5);

        onering(1) = f;
        e  = mesh.edges(f * 3);
        for v = 1:v1 - 2
            e = e.next.pair;
            onering(1 + v) = findedge(e);
        end
        e  = mesh.edges(f * 3).next;
        for v = 1:v2 - 2
            e = e.next.pair;
            onering(v1 - 1 + v) = findedge(e);
        end
        e  = mesh.edges(f * 3).next.next;
        for v = 1:v3 - 2
            e = e.next.pair;
            onering(v1 + v2 - 3 + v) = findedge(e);
        end
        
        for c = 1:ptsperface
            cd = P + (f - 1) * ptsperface + c;
            [cand_limits(cd, :) deriv] = ...
                mesh.limitevaluation(mesh.edges(cand_ixs(cd)), ...
                    cand_uvs(cd, 1), cand_uvs(cd, 2));
            cand_derivs(cd, :, 1) = deriv(2, :);
            cand_derivs(cd, :, 2) = [ -1/sqrt(3) 0 1/sqrt(3) ] * deriv;
        end
                        
        for o = 1:length(onering)
            if ~faces_done(f, onering(o))
                for c = 1:ptsperface
                    for d = 1:ptsperface
                        cd  = P + (f - 1) * ptsperface + c;
                        cd2 = P + (onering(o) - 1) * ptsperface + d;
                        cand_dists(cd, cd2)  = ...
                            mesh.distbetween(mesh.edges(cand_ixs(cd)), ...
                                             cand_uvs(cd, 1), ...
                                             cand_uvs(cd, 2), ...
                                             mesh.edges(cand_ixs(cd2)), ...
                                             cand_uvs(cd2, 1), ...
                                             cand_uvs(cd2, 2));
                        cand_dists(cd2, cd)  = cand_dists(cd, cd2);
                    end
                end
                faces_done(f, onering(o)) = true;
                faces_done(onering(o), f) = true;
            end
        end
    end
    fprintf(1, '\n\n');
    
    function facenum = findedge(e)
        if any(e.index_in_mesh == face_edge_table)
            facenum = find(e.index_in_mesh == face_edge_table);
        elseif any(e.next.index_in_mesh == face_edge_table)
            facenum = find(e.next.index_in_mesh == face_edge_table);
        else
            facenum = find(e.next.next.index_in_mesh == face_edge_table);
        end
    end
end


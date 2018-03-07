function [ vertices faces ] = trimesh_fromply(plyfile)

% TRIMESH_FROMPLY  Load vertices and faces from .PLY file

    scalar_properties(1).element  = 'vertex';
    scalar_properties(1).property = 'x';
    scalar_properties(1).callback = @(x) addvertex(x, 1);
    scalar_properties(2).element  = 'vertex';
    scalar_properties(2).property = 'y';
    scalar_properties(2).callback = @(x) addvertex(x, 2);
    scalar_properties(3).element = 'vertex';
    scalar_properties(3).property = 'z';
    scalar_properties(3).callback = @(x) addvertex(x, 3);

    list_properties(1).element = 'face';
    list_properties(1).property = 'vertex_indices';
    list_properties(1).start = @(x) resetfaceposition();
    list_properties(1).item  = @addvindex;

    vertices     = [];
    faces        = [];
    faceposition = 1;

    parse_ply(plyfile, scalar_properties, list_properties);

    faces = faces + 1;
    
    function addvertex(val, ix)
        if ix == 1
            vertices(end + 1, ix) = val;
        else
            vertices(end, ix) = val;
        end
    end

    function addvindex(val)
        if faceposition == 1
            faces(end + 1, faceposition) = val;
        else
            faces(end, faceposition) = val;
        end
        faceposition = faceposition + 1;
    end

    function resetfaceposition
        faceposition = 1;
    end
end
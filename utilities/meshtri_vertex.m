classdef meshtri_vertex < handle
    % MESHTRI_VERTEX  Class used to represent the vertex of a triangle mesh.
    %
    %   See also MESHTRI
    
    properties
        edge
        definition
        valency
        index_in_mesh
        
        newVert
    end
    
    methods
        function obj = meshtri_vertex(totalverts, vertexnum)
            if nargin > 0
                obj.definition = zeros(1, totalverts);
            end
            
            if nargin > 1
                obj.definition(vertexnum) = 1;
                obj.index_in_mesh         = vertexnum;
            end
        end
        
        function setValency(obj)
            obj.valency = 1;
            iterEdge = obj.edge.pair.next;
            
            while iterEdge ~= obj.edge
                obj.valency = obj.valency + 1;
                iterEdge = iterEdge.pair.next;
            end
        end
        
        function indices = adjVerts(obj)
            indices    = zeros(1, obj.valency);
            indices(1) = obj.edge.vert.index_in_mesh;
            iterEdge   = obj.edge.pair.next;
            ind = 1;
            
            while iterEdge ~= obj.edge
                ind = ind + 1;
                indices(ind) = iterEdge.vert.index_in_mesh;
                iterEdge = iterEdge.pair.next;
            end
        end
        
        function setNewVert(obj)
            obj.newVert = meshtri_vertex(length(obj.definition));
            
            alpha = (0.375 + cos(2 * pi / obj.valency) / 4) ^ 2 + 0.375;
            delta = (1 - alpha) / obj.valency;
            
            done     = false;
            iterEdge = obj.edge;
            while ~done
                obj.newVert.definition = obj.newVert.definition + ...
                    iterEdge.vert.definition * delta;
                iterEdge = iterEdge.pair.next;
                done = (iterEdge == obj.edge);
            end
            
            obj.newVert.definition = obj.newVert.definition + ...
                obj.definition * alpha;
            
            obj.newVert.definition = obj.newVert.definition / ...
                sum(obj.newVert.definition);
        end
    end

end


classdef meshtri_halfedge < handle
    % MESHTRI_HALFEDGE  Class used to represent a half-edge
    %   in a triangle mesh.
    %
    %   See also MESHTRI
    
    properties
        pair
        next
        vert
        edge
        index_in_mesh
    end
    
    methods
        function obj = meshtri_halfedge(vert)
            if nargin > 0
                obj.vert = vert;
            end
        end
        
        function setNewVert(obj)
            if isempty(obj.edge.newVert)
                obj.edge.newVert = ...
                    meshtri_vertex(length(obj.vert.definition));
                
                newVert = obj.edge.newVert;
                newVert.definition = newVert.definition + ...
                    0.375 * obj.vert.definition;
                newVert.definition = newVert.definition + ...
                    0.375 * obj.pair.vert.definition;
                
                newVert.definition = newVert.definition + ...
                    0.125 * obj.next.vert.definition;
                newVert.definition = newVert.definition + ...
                    0.125 * obj.pair.next.vert.definition;
                
                newVert.definition = newVert.definition / ...
                                     sum(newVert.definition);
            end
        end
        
        function res = isInLoop(obj, edge)
            res = edge == obj || edge == obj.next || edge == obj.next.next;
        end
    end
    
end


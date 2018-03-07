function parse_ply(plyfile, ...
                   scalar_prop_callbacks, ...
                   list_prop_callbacks, ...
                   element_callbacks)

% PARSE_PLY  Simple parser for .PLY files in ASCII format
%
%   See also TRIMESH_FROMPLY
               
%#ok<*STTOK>
%#ok<*AGROW>

if nargin < 2, scalar_prop_callbacks = []; end
if nargin < 3, list_prop_callbacks = [];   end
if nargin < 4, element_callbacks = [];     end

fid = fopen(plyfile);

magic = fgetl(fid);
if ~strcmp(magic, 'ply')
    error('Magic bytes are incorrect for PLY file');
end

format = fgetl(fid);
if ~strcmp(format, 'format ascii 1.0');
    error('Format unsupported');
end

file_definition = [];
while ~feof(fid)
    line = fgetl(fid);
    [tok, line] = strtok(line);
    
    if strcmp(tok, 'comment')
        continue;
    elseif strcmp(tok, 'element')
        [el_name, line] = strtok(line);
        el_count        = strtok(line);
        
        file_definition(end + 1).name   = el_name;
        file_definition(end).count      = str2double(el_count);
        file_definition(end).properties = [];
    elseif strcmp(tok, 'property')
        [tok, line] = strtok(line);
        
        file_definition(end).properties(end + 1).type = tok;
        if strcmp(tok, 'list')
            [count_type, line] = strtok(line);
            [item_type,  line] = strtok(line);
            
            file_definition(end).properties(end).count_type = count_type;
            file_definition(end).properties(end).item_type  = item_type;
        end
        prop_name = strtok(line);
        file_definition(end).properties(end).name = prop_name;
    elseif strcmp(tok, 'end_header')
        break;
    else
        error('Unknown keyword');
    end
end

for element = 1:length(file_definition)
    for instance = 1:file_definition(element).count
        % Start element callbacks
        for i = 1:length(element_callbacks)
            if strcmp(element_callbacks(i).element, ...
                      file_definition(element).name)
                if isfield(element_callbacks(i), 'start')
                    element_callbacks(i).start();
                end
            end
        end
        
        element_data = fgetl(fid);
        for i = 1:length(file_definition(element).properties)
            if strcmp(file_definition(element).properties(i).type, 'list')
                [count, element_data] = strtok(element_data); 
                count = str2double(count);
                
                % Start list callbacks
                for j = 1:length(list_prop_callbacks)
                    if strcmp(list_prop_callbacks(j).element,   ...
                              file_definition(element).name) && ...
                       strcmp(list_prop_callbacks(j).property,  ...
                              file_definition(element).properties(i).name)
                        if isfield(list_prop_callbacks(j), 'start')
                            list_prop_callbacks(j).start(count);
                        end
                    end
                end
                
                for item_num = 1:count
                    [item, element_data] = strtok(element_data); 
                    item = str2double(item);
                    
                    % List item callbacks
                    for j = 1:length(list_prop_callbacks)
                        if strcmp(list_prop_callbacks(j).element,   ...
                                  file_definition(element).name) && ...
                           strcmp(list_prop_callbacks(j).property,  ...
                                  file_definition(element).properties(i).name)
                            if isfield(list_prop_callbacks(j), 'item')
                                list_prop_callbacks(j).item(item);
                            end
                        end
                    end
                end
                
                % End list callbacks
                for j = 1:length(list_prop_callbacks)
                    if strcmp(list_prop_callbacks(j).element,   ...
                              file_definition(element).name) && ...
                       strcmp(list_prop_callbacks(j).property,  ...
                              file_definition(element).properties(i).name)
                        if isfield(list_prop_callbacks(j), 'end')
                            list_prop_callbacks(j).end();
                        end
                    end
                end
            else
                [property, element_data] = strtok(element_data); 
                property = str2double(property);
                
                % Scalar property callbacks
                for j = 1:length(scalar_prop_callbacks)
                    if strcmp(scalar_prop_callbacks(j).element,  ...
                              file_definition(element).name) &&  ...
                       strcmp(scalar_prop_callbacks(j).property, ...
                              file_definition(element).properties(i).name)
                        if isfield(scalar_prop_callbacks(j), 'callback')
                            scalar_prop_callbacks(j).callback(property);
                        end
                    end
                end
            end
        end
        
        % End element callbacks
        for i = 1:length(element_callbacks)
            if strcmp(element_callbacks(i).element, ...
                      file_definition(element).name)
                if isfield(element_callbacks(i), 'end')
                    element_callbacks(i).end();
                end
            end
        end
    end
end

fclose(fid);

end


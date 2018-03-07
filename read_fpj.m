function [ plyfile images ] = read_fpj(file)

% READ_FPJ  Read Forms project file.
%
%   Call through READ_PROJECT to set up a project struct used in FORMS.
%
%   See also READ_PROJECT, FORMS

fid = fopen(file);

magic = fgetl(fid);
assert(strcmp(magic, 'fpj'));
format = fgetl(fid);
assert(strcmp(format, 'format ascii 1.0'));

rootdir = fileparts(file);
plyfile = fullfile(rootdir, fgetl(fid));

images = [];
while ~feof(fid)
    image.filename         = fullfile(rootdir, fgetl(fid));
    
    % Size, silhouette and constraints are actually given in measure units
    % (1/96th of an inch) rather than pixels, but that's irrelevant: we're
    % about to normalize anyway
    sizestring             = fgetl(fid);
    imagesize              = textscan(sizestring, '%f', 'delimiter', ',');
    imagesize              = imagesize{1};
    
    image.pathstring       = fgetl(fid);
    [image.points ...
     image.normalsLeft]    = parse_fpjpath(image.pathstring);
    
    normalizescale         = 2 / imagesize(1);
    image.points(:, 1, :)  = (normalizescale * image.points(:, 1, :)) - 1;
    image.points(:, 2, :)  = (normalizescale * imagesize(2) * 0.5) - ...
                             (normalizescale * image.points(:, 2, :));
    
    image.constraintstring = fgetl(fid);
    [p3d p2d onsil]        = parse_fpjconstraints(image.constraintstring);
    
    p2d(:, 1)              = (normalizescale * p2d(:, 1)) - 1;
    p2d(:, 2)              = (normalizescale * imagesize(2) * 0.5) - ...
                             (normalizescale * p2d(:, 2));
    
    image.constraints3d    = p3d;
    image.constraints2d    = p2d;
    image.constraintsonsil = onsil;
    
    image.transformstring  = fgetl(fid);
    trans_vals = textscan(image.transformstring, '%f', 'delimiter', ',');
    trans_vals = trans_vals{1};
    image.translate        = eye(4);
    image.translate(1, 4)  = trans_vals(4);
    image.translate(2, 4)  = trans_vals(5);
    image.rotate           = ...
        [ expm([    0 -trans_vals(3)  trans_vals(2) ;
                 trans_vals(3)     0 -trans_vals(1) ;
                -trans_vals(2)  trans_vals(1)     0 ]) zeros(3, 1) ;
                                                       zeros(1, 3) 1 ];
    image.scale            = trans_vals(6) * eye(4);
    image.scale(4, 4)      = 1;
    image.transform        = image.translate * image.rotate * image.scale;
    
    images = [ images image ];                                 %#ok<AGROW>
end

fclose(fid);

end
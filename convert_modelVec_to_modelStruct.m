function model = convert_modelVec_to_modelStruct(modelVec, project)

% Function for converting the vector version of the model (modelVec, the
% output of Cashman & Fitzgibbon's original version of forms.m), to the
% structure version (as in the current version of forms.m)

% Jason Manley, Nov 2017

model = struct();

model.n = modelVec(1);
model.M = modelVec(2);
model.S = modelVec(3:2+model.n);
model.P = length(project.mesh.verts);
model.shapemodes = reshape(modelVec(end - 3 * model.P * (model.M + 1) + 1:end), ...
    3 * model.P, model.M + 1);

model.shapevars = cell(model.n,1);
model.scale     = zeros(model.n,1);
model.rotate    = cell(model.n,1);
model.translate = cell(model.n,1);

for i=1:model.n
    vs         = 5 * sum(model.S(1:i - 1)) + (i - 1) * (7 + model.M) + model.n + 2;
    tvs        = vs + 5 * model.S(i) + 1;
    tve        = vs + 5 * model.S(i) + 7;
    mvs        = vs + 5 * model.S(i) + 8;
    mve        = vs + 5 * model.S(i) + 7 + model.M;
    v          = modelVec(tvs:tvs + 2);
    rotate   = [expm([    0 -v(3)  v(2) ;
        v(3)     0 -v(1) ;
        -v(2)  v(1)     0 ]) zeros(3, 1) ;
        zeros(1, 3) 1];
    translate  = [ 1 0 0 modelVec(tvs + 3) ;
        0 1 0 modelVec(tvs + 4) ;
        0 0 1 modelVec(tvs + 5) ;
        0 0 0 1              ];
    
    model.shapevars{i} = modelVec(mvs:mve);
    model.scale(i)     = modelVec(tve);
    model.rotate{i}    = rotate;
    model.translate{i} = translate;
end


end
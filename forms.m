function [model, modelVec] = forms(project, parameters)

% FORMS  Finds morphable model from silhouettes given in 'project'.
%
%   Create a 'project' struct using READ_PROJECT for original
%   Cashman/Fitzgibbon datasets.
%
%   See README for a description of curating a brand new 'project' struct.
%
%   parameters.M           is the number of deformation basis shapes (D in the paper)
%   parameters.xi_0        is the weight on E^tp_0
%   parameters.xi_def      is the weight on E^tp_m for m > 0
%   parameters.sigma_norm  is the noise estimate on normals
%   parameters.sigma_con   is the noise estimate on point constraints
%   parameters.sigma_sil   is the noise estimate on silhouettes
%
%   See also READ_PROJECT, ALIGN_PROJECT

% --
% -- Initialization
% --

% updated Jason Manley, Dec 2017
% Questions? jmanley@rockefeller.edu

if nargin < 2
    parameters = makeModelParameters;
end

% User-defined connectivity (mesh)
mesh             = project.mesh;
% Thin plate energy
thinplate        = mesh.thinplate();
tplatesqrt       = real(sqrtm(thinplate / 2));
tplatesqrt_3     = blkdiag(tplatesqrt, tplatesqrt, tplatesqrt);
% Set this to false for slow, checked Jacobian computation
jacobmult_on     = true;

% Number of points in the mesh
P                = length(mesh.verts);
% Number of images
n                = length(project.images);

% Weight on shape coefficient regularisation terms
beta             = 0.5;
% Weight on contour generator continuity terms
gamma            = 1 / 128;

% Number of constraints in each image
K                = zeros(n, 1);
% Weight on each constraint point
conwt            = cell(n, 1);
% Resolution of silhouette in each image
S                = zeros(n, 1);
% Number of free contour generator points
T                = zeros(n, 1);
for i = 1:n
    K(i)         = length(project.images(i).constraints3d);
    conwt{i}     = ones(K(i), 1) / parameters.sigma_con;
    S(i)         = 125;
    T(i)         = S(i) - sum(project.images(i).constraintsonsil);
end

% --
% -- A large pile of index vectors to make life easier
% --

% varstart(i) is 1 less than the first variable for image i
vs               = zeros(n + 1, 1);
% silvarstart(i):silvarend(i) are the silhouette vars for image i
svs              = zeros(n    , 1);
sve              = zeros(n    , 1);
% rotvarstart(i):rotvarend(i) are the rotation vars for image i
rvs              = zeros(n    , 1);
rve              = zeros(n    , 1);
% transvarstart(i):transvarend(i) are the translation variables for
% image i
tvs              = zeros(n    , 1);
tve              = zeros(n    , 1);
% scalevar(i) is the scale variable for image i
sv               = zeros(n    , 1);
% modevarstart(i):modevarend(i) are the shape variables for image i
mvs              = zeros(n    , 1);
mve              = zeros(n    , 1);

% eqstart(i) is 1 less than the first equation for image i
es               = zeros(n + 1, 1);
% sileqstart(i):sileqend(i) are the silhouette equations for image i
ses              = zeros(n    , 1);
see              = zeros(n    , 1);
% coneqstart(i):coneqend(i) are the silhouette continuity equations for
% image i
ces              = zeros(n    , 1);
cee              = zeros(n    , 1);
% normeqstart(i):normeqend(i) are the silhouette normal equations
% for image i
nes              = zeros(n    , 1);
nee              = zeros(n    , 1);
% constxeqstart(i):constxeqend(i) are the constraint equations for
% the x coordinates on image i
cxes             = zeros(n    , 1);
cxee             = zeros(n    , 1);
% constyeqstart(i):constyeqend(i) are the constraint equations for
% the y coordinates on image i
cyes             = zeros(n    , 1);
cyee             = zeros(n    , 1);
% modeeqstart(i):modeeqend(i) are the shape coefficient regularisation
% equations for image i
mes              = zeros(n    , 1);
mee              = zeros(n    , 1);

% --
% -- Basis vectors and globals for update_sil_surface
% --

u_basis          = complex(1, 0);
v_basis          = exp(complex(0, pi / 3));
change_basis     = [ real(u_basis) real(v_basis) ;
    imag(u_basis) imag(v_basis) ];
vertex_radius    = 0.49; % (must be strictly less than 0.5)
circ_bot         = vertex_radius * u_basis;
circ_top         = vertex_radius * v_basis;

% --
% -- Constraint points
% --

constAeq         = cell(n, 1);
constbeq         = cell(n, 1);
for i = 1:n
    constverts      = project.images(i).constraints3d;
    constAeq{i}     = zeros(2 * K(i), 2 * P);
    
    for k = 1:K(i)
        e           = mesh.verts(constverts(k)).edge.next;
        limitpoint  = mesh.limitevaluation(e, 0, 0);
        
        % X
        constAeq{i}(       k,     1:    P) = limitpoint;
        % Y
        constAeq{i}(K(i) + k, P + 1:2 * P) = limitpoint;
    end
    constbeq{i}     = project.images(i).constraints2d(:);
end

% --
% -- Initialization for the optimizer
% --

if parameters.M> 0, sM = 1; else sM = 0; end
problem.x0          = zeros(2 * sum(T) + (7 + sM) * n, 1);
last_sil_uvs        = cell(n, 1);
for i = 1:n
    vs(i + 1)       = vs(i) + 2 * T(i) + 7 + sM;
    last_sil_uvs{i} = zeros(2 * T(i), 1);
    problem.x0(vs(i) + 1:vs(i + 1) - sM) ...
        = [ last_sil_uvs{i} ; ...   Silhouette variables
        zeros(6, 1) ; 1 ];  %   Transform variables
end
problem.x0          = [ problem.x0 ; project.vertices(:) ];
if sM > 0
    problem.x0      = [ problem.x0 ; zeros(3 * P, 1) ];
end
silJ                = cell(n, 1);
varsJ               = cell(n, 1);
meshJ               = cell(n, 1);

% --
% -- Silhouette points
% --

sil_pts                  = cell(n, 1);
sil_normals              = cell(n, 1);
fixed_sil_pts            = cell(n, 1);
fixed_indices            = cell(n, 1);
mu                       = cell(n, 1);
sil_triangles            = cell(n, 1);
sil_barycentric          = cell(n, 1);
last_sil_barycen         = cell(n, 1);
last_sil_tris            = cell(n, 1);
for i = 1:n
    fixed_sil_pts{i}     = false(S(i), 1);
    silhouette           = project.images(i).points;
    sil_params           = sil_sample(S(i), silhouette);
    sil_pts{i}           = zeros(S(i), 2);
    sil_normals{i}       = zeros(S(i), 2);
    
    for s = 1:S(i)
        zerotangent = true;
        while zerotangent
            seg              = floor(sil_params(s));
            sil_pts{i}(s, :) = sil_evalbezier(...
                silhouette(:, :, 1 + seg), ...
                sil_params(s) - seg);
            tan_pts          = 3 * (silhouette(2:end, :, 1 + seg) - ...
                silhouette(1:end - 1, :, 1 + seg));
            tangent          = sil_evalbezier(tan_pts, ...
                sil_params(s) - seg);
            
            zerotangent      = (norm(tangent, 2) == 0);
            % Next place to try sampling the silhouette, if it turned
            % out that the current place has zero first derivative.
            % (sil_params is not used again, so it's safe to modify
            % it).
            sil_params(s)    = sil_params(s) + 1e-4;
        end
        sil_tan = tangent / norm(tangent, 2);
        
        if project.images(i).normalsLeft
            sil_normals{i}(s, :) = [ -sil_tan(2)  sil_tan(1) ];
        else
            sil_normals{i}(s, :) = [  sil_tan(2) -sil_tan(1) ];
        end
    end
    
    mu{i}            = zeros(sum(project.images(i).constraintsonsil), 1);
    fixed_indices{i} = zeros(sum(project.images(i).constraintsonsil), 1);
    f  = 0;
    for k = 1:K(i)
        if project.images(i).constraintsonsil(k)
            % Snap constraint points on silhouette to nearest
            % silhouette point
            closest = sil_pts{i} - ...
                repmat(project.images(i).constraints2d(k, :), ...
                S(i), 1);
            [~, ind]              = min(sum(closest .^ 2, 2));
            constbeq{i}(k)        = sil_pts{i}(ind, 1);
            constbeq{i}(K(i) + k) = sil_pts{i}(ind, 2);
            assert(fixed_sil_pts{i}(ind) == false);
            fixed_sil_pts{i}(ind) = true;
            f                     = f + 1;
            mu{i}(f)              = project.images(i).constraints3d(k);
            fixed_indices{i}(f)   = ind;
        end
    end
    
    [~, correct_mu] = sort(fixed_indices{i});
    mu{i} = mu{i}(correct_mu);
end

init_sigma_sil   = parameters.sigma_sil;
init_sigma_norm  = parameters.sigma_norm / 3;
init_gamma       = gamma * 9;
for cM = sM:parameters.M
    fprintf(1, '\nFitting dimension %d out of %d...\n', cM, parameters.M);
    for i = 1:n
        vs(i + 1) = vs(i) + 2 * T(i) + 7 + cM;
        svs(i)    = vs(i) + 1;
        sve(i)    = vs(i) + 2 * T(i);
        rvs(i)    = vs(i) + 2 * T(i) + 1;
        rve(i)    = vs(i) + 2 * T(i) + 3;
        tvs(i)    = vs(i) + 2 * T(i) + 4;
        tve(i)    = vs(i) + 2 * T(i) + 6;
        sv(i)    = vs(i) + 2 * T(i) + 7;
        mvs(i)    = vs(i) + 2 * T(i) + 8;
        mve(i)    = vs(i) + 2 * T(i) + 7 + cM;
        
        es(i + 1) = es(i) + 5 * T(i) + S(i) + 2 * K(i) + cM;
        ses(i)    = es(i) + 1;
        see(i)    = es(i) + 2 * T(i);
        ces(i)    = es(i) + 2 * T(i) + 1;
        cee(i)    = es(i) + 2 * T(i) + S(i);
        nes(i)    = es(i) + 2 * T(i) + S(i) + 1;
        nee(i)    = es(i) + 5 * T(i) + S(i);
        cxes(i)   = es(i) + 5 * T(i) + S(i) + 1;
        cxee(i)   = es(i) + 5 * T(i) + S(i) +     K(i);
        cyes(i)   = es(i) + 5 * T(i) + S(i) +     K(i) + 1;
        cyee(i)   = es(i) + 5 * T(i) + S(i) + 2 * K(i);
        mes(i)    = es(i) + 5 * T(i) + S(i) + 2 * K(i) + 1;
        mee(i)    = es(i) + 5 * T(i) + S(i) + 2 * K(i) + cM;
    end
    
    fprintf(1, '\nFinding %d contour generators\n', n);
    modes = reshape(problem.x0(vs(n + 1) + 1:end), 3 * P, cM + 1);
    for i = 1:n
        fprintf(1, '%d ', i);
        v                = problem.x0(rvs(i):rve(i));
        rotscale         = [ problem.x0(sv(i)) .* ...
            expm([    0 -v(3)  v(2) ;
            v(3)     0 -v(1) ;
            -v(2)  v(1)     0 ]) zeros(3, 1) ;
            zeros(1, 3) 1];
        translate        = [ 1 0 0 problem.x0(tvs(i)    ) ;
            0 1 0 problem.x0(tvs(i) + 1) ;
            0 0 1 problem.x0(tvs(i) + 2) ;
            0 0 0 1                ];
        DT               = translate * ...
            project.images(i).transform * rotscale;
        pointvars        = reshape(modes * ...
            [1 ; problem.x0(mvs(i):mve(i))], ...
            P, 3);
        transformed      = [ pointvars ones(P, 1) ] * DT';
        transformed      = transformed(:, 1:3);
        cand_norms       = cross(project.cand_derivs(:, :, 2) * ...
            transformed, ...
            project.cand_derivs(:, :, 1) * ...
            transformed);
        cand_norms       = cand_norms ./ ...
            repmat(sqrt(sum(cand_norms .^ 2, 2)), 1, 3);
        transformed      = transformed(:, 1:2);
        
        if any(fixed_sil_pts{i})
            sil_selcands = zeros(1, S(i));
            nf           = sum(fixed_sil_pts{i});
            fixed_cands  = find(fixed_sil_pts{i});
            for f = 1:nf
                if f == nf, g = 1; else g = f + 1; end
                
                lin_path = sil_conspreimage(sil_pts{i}, ...
                    sil_normals{i}, init_sigma_sil, ...
                    init_gamma, init_sigma_norm, ...
                    project.cand_limits * transformed, ...
                    project.cand_dists, cand_norms, ...
                    [mu{i}(f) mu{i}(f)], [mu{i}(g) mu{i}(g)], ...
                    [fixed_cands(f) fixed_cands(g)], false);
                
                if f < nf
                    sil_selcands(fixed_cands(f):fixed_cands(g)) = ...
                        lin_path;
                else
                    sil_selcands(fixed_cands(f):end) = ...
                        lin_path(1:S(i) - fixed_cands(f) + 1);
                    sil_selcands(1:fixed_cands(g)) = ...
                        lin_path(S(i) - fixed_cands(f) + 2:end);
                end
            end
        else
            sil_selcands = sil_circpreimage(sil_pts{i}, ...
                sil_normals{i}, init_sigma_sil, ...
                init_gamma, init_sigma_norm, ...
                project.cand_limits * transformed, ...
                project.cand_dists, cand_norms);
        end
        sil_triangles{i}     = project.cand_ixs(sil_selcands);
        sil_barycentric{i}   = project.cand_uvs(sil_selcands, :);
        
        for s = 1:S(i)
            u = sil_barycentric{i}(s, 1);
            v = sil_barycentric{i}(s, 2);
            w = 1 - sum(sil_barycentric{i}(s, :));
            [~, tripos] = max([u w v]);
            
            % This function assumes that 'tripos' is 2. If it isn't,
            % rotate round.
            if tripos == 1
                sil_barycentric{i}(s, 1) = sil_barycentric{i}(s, 2);
                sil_barycentric{i}(s, 2) = w;
                sil_triangles{i}(s) = ...
                    mesh.edges(sil_triangles{i}(s)).next.next.index_in_mesh;
            elseif tripos == 3
                sil_barycentric{i}(s, 2) = sil_barycentric{i}(s, 1);
                sil_barycentric{i}(s, 1) = w;
                sil_triangles{i}(s) = ...
                    mesh.edges(sil_triangles{i}(s)).next.index_in_mesh;
            end
        end
        
        % Pull all non-fixed points slightly away from vertices
        ptsatverts = sum(abs(sil_barycentric{i}), 2) == 0 & ...
            ~fixed_sil_pts{i};
        sil_barycentric{i}(ptsatverts, :) = 1e-2;
        last_sil_barycen{i}  = sil_barycentric{i};
        last_sil_tris{i}     = sil_triangles{i};
    end
    fprintf(1, '\n');
    
    if cM > 0
        for i = 1:n
            problem.x0(mvs(i) + cM - 1) = 1;
        end
    end
    
    % --
    % -- Invoke lsqnonlin optimization
    % --
    problem.objective   = @calc_energy;
    problem.options     = optimset( ...
        'Jacobian',         'on'                 ...
        , 'PreCondBandwidth', Inf                  ...
        , 'Diagnostics',      'on'                 ...
        , 'Display',          'iter'               ...
        , 'OutputFcn',        @newiteration        ...
        , 'TolFun',           1e-4                 ...
        , 'MaxIter',          100                  ...
        );
    if jacobmult_on
        problem.options = optimset(problem.options ...
            , 'JacobMult',        @jacobmult           ...
            );
    else
        problem.options = optimset(problem.options ...
            , 'DerivativeCheck', 'on'                  ...
            );
    end
    problem.solver                            = 'lsqnonlin';
    [result,resnorm,residual,exitflag,output] = lsqnonlin(problem);
    
    model = zeros(5 * sum(S) + n * (7 + cM), 1);
    currpos = 1;
    for i = 1:n
        silhouette_vars = [ sil_triangles{i} sil_barycentric{i} ]';
        model(currpos:currpos + 3 * S(i) - 1) = silhouette_vars(:);
        model(currpos + 3 * S(i):currpos + 5 * S(i) - 1) = sil_pts{i}';
        model(currpos + 5 * S(i):currpos + 5 * S(i) + 6 + cM) = ...
            result(rvs(i):mve(i));
        currpos = currpos + 5 * S(i) + 7 + cM;
    end
    model = [ n ; cM ; S ; model ; result(vs(n + 1) + 1:end) ];         %#ok<AGROW>
    
    % --
    % -- Initialization to learn the next shape mode
    % --
    
    problem.x0          = zeros(vs(n + 1) + n + (cM + 2) * 3 * P, 1);
    last_sil_uvs        = cell(n, 1);
    for i = 1:n
        last_sil_uvs{i} = zeros(2 * T(i), 1);
        problem.x0(svs(i) + i - 1:sve(i) + i - 1) ...
            = last_sil_uvs{i};
        problem.x0(rvs(i) + i - 1:mve(i) + i - 1) ...
            = result(rvs(i):mve(i));
    end
    problem.x0(vs(n + 1) + n + 1:vs(n + 1) + n + (cM + 1) * 3 * P) ...
        = result(vs(n + 1) + 1:end);
    
    % Make the initialization problem converge to the
    % geometry-modification problem as more shape modes are added
    init_sigma_sil  = 2 / (1 / parameters.sigma_sil  + 1 / init_sigma_sil );
    init_sigma_norm = 2 / (1 / parameters.sigma_norm + 1 / init_sigma_norm);
    init_gamma      = (((sqrt(2*gamma) + sqrt(2*init_gamma))/2) ^ 2)/2;
end

% Convert vector model into a more sensible output

modelVec = model;
model = struct();

model.n = n;
model.M = parameters.M;
model.parameters = parameters;
model.parameters.beta = beta;
model.parameters.gamma = gamma;
model.S = S;
model.T = T;
model.K = K;
model.P = P;
model.conwt = conwt;
model.shapemodes = reshape(modelVec(end - 3 * P * (model.M + 1) + 1:end), ...
    3 * P, model.M + 1);
model.resnorm   = resnorm;
model.residual  = residual;
model.exitflag  = exitflag;
model.output    = output;

model.shapevars = cell(n,1);
model.scale     = zeros(n,1);
model.rotate    = cell(n,1);
model.translate = cell(n,1);

for i=1:n
    vs         = 5 * sum(S(1:i - 1)) + (i - 1) * (7 + model.M) + n + 2;
    tvs        = vs + 5 * S(i) + 1;
    tve        = vs + 5 * S(i) + 7;
    mvs        = vs + 5 * S(i) + 8;
    mve        = vs + 5 * S(i) + 7 + model.M;
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


    function stop = newiteration(vars, ~, state)
        if strcmp(state, 'iter')
            for i = 1:n
                sil_barycentric{i}  = last_sil_barycen{i};
                sil_triangles{i}    = last_sil_tris{i};
                update_sil_surface(i, vars(svs(i):sve(i)));
                last_sil_barycen{i} = sil_barycentric{i};
                last_sil_tris{i}    = sil_triangles{i};
                last_sil_uvs{i}     = vars(svs(i):sve(i));
            end
        end
        stop = false;
    end

    function W = jacobmult(Jdata, X, flag)
        vars = Jdata(end, 1:2 * sum(T) + n * (7 + cM) + 3 * P * (cM + 1));
        
        meansc = 0;
        for i = 1:n
            bs       = (i - 1) * (9 + cM + 3 * P) + 1;
            be       =       i * (9 + cM + 3 * P);
            
            silJ{i}  = Jdata(1:2 * S(i) + 5 * T(i), bs:bs + 1);
            varsJ{i} = Jdata(1:    S(i) + 5 * T(i) + 2 * K(i), ...
                bs + 2:bs + 8 + cM);
            meshJ{i} = Jdata(1:    S(i) + 5 * T(i) + 2 * K(i), ...
                bs + 9 + cM:be);
            
            meansc   = meansc + vars(sv(i));
        end
        meansc = meansc / n;
        
        if flag == 0
            flag = 1;
            both = 1;
        else
            both = 0;
        end
        
        if flag > 0
            Y = zeros(es(n + 1) + 3 * P * (cM + 1), size(X, 2));
            
            for i = 1:n
                esi = es(i) + 1;
                esn = es(i + 1) - cM;
                Y(esi:esn, :) = varsJ{i} * X(rvs(i):mve(i), :);
                Y(esi:esn, :) = Y(esi:esn, :) + ...
                    meshJ{i} * X(vs(n + 1) + 1:vs(n + 1) + 3 * P, :);
                for m = 1:cM
                    Y(esi:esn, :) = Y(esi:esn, :) + meshJ{i} * ...
                        vars(mvs(i) + m - 1) * X(vs(n + 1) + 1 + ...
                        m * 3 * P:vs(n + 1) + (m + 1) * 3 * P, :);
                end
                
                Y(ses(i):see(i) - T(i), :) = Y(ses(i):see(i) - T(i), :) ...
                    + diag(silJ{i}(1:T(i), 1)) * X(svs(i):sve(i) - T(i), :);
                Y(ses(i):see(i) - T(i), :) = Y(ses(i):see(i) - T(i), :) ...
                    + diag(silJ{i}(1:T(i), 2)) * X(svs(i) + T(i):sve(i), :);
                Y(ses(i) + T(i):see(i), :) = Y(ses(i) + T(i):see(i), :) ...
                    + diag(silJ{i}(T(i) + 1:2 * T(i), 1)) ...
                    * X(svs(i):sve(i) - T(i), :);
                Y(ses(i) + T(i):see(i), :) = Y(ses(i) + T(i):see(i), :)...
                    + diag(silJ{i}(T(i) + 1:2 * T(i), 2)) ...
                    * X(svs(i) + T(i):sve(i), :);
                
                D = diag(silJ{i}(2 * T(i) + 1:2 * T(i) + S(i), 1)) + ...
                    diag(silJ{i}(2 * T(i) + S(i) + 1:...
                    2 * T(i) + 2 * S(i) - 1, 1), 1);
                D(S(i), 1) = silJ{i}(2 * T(i) + 2 * S(i), 1);
                D(:, fixed_sil_pts{i}) = [];
                Y(ces(i):cee(i), :) = Y(ces(i):cee(i), :) ...
                    +  D * X(svs(i):sve(i) - T(i), :);
                D = diag(silJ{i}(2 * T(i) + 1:2 * T(i) + S(i), 2)) + ...
                    diag(silJ{i}(2 * T(i) + S(i) + 1:...
                    2 * T(i) + 2 * S(i) - 1, 2), 1);
                D(S(i), 1) = silJ{i}(2 * T(i) + 2 * S(i), 2);
                D(:, fixed_sil_pts{i}) = [];
                Y(ces(i):cee(i), :) = Y(ces(i):cee(i), :) ...
                    + D * X(svs(i) + T(i):sve(i), :);
                
                Y(nes(i):nee(i) - 2 * T(i), :) = Y(nes(i):nee(i) - 2 * T(i), :) ...
                    + diag(silJ{i}(2 * T(i) + 2 * S(i) + 1:3 * T(i) + 2 * S(i), 1)) ...
                    * X(svs(i):sve(i) - T(i), :);
                Y(nes(i):nee(i) - 2 * T(i), :) = Y(nes(i):nee(i) - 2 * T(i), :) ...
                    + diag(silJ{i}(2 * T(i) + 2 * S(i) + 1:3 * T(i) + 2 * S(i), 2)) ...
                    * X(svs(i) + T(i):sve(i), :);
                Y(nes(i) + T(i):nee(i) - T(i), :) = Y(nes(i) + T(i):nee(i) - T(i), :) ...
                    + diag(silJ{i}(3 * T(i) + 2 * S(i) + 1:4 * T(i) + 2 * S(i), 1)) ...
                    * X(svs(i):sve(i) - T(i), :);
                Y(nes(i) + T(i):nee(i) - T(i), :) = Y(nes(i) + T(i):nee(i) - T(i), :)...
                    + diag(silJ{i}(3 * T(i) + 2 * S(i) + 1:4 * T(i) + 2 * S(i), 2)) ...
                    * X(svs(i) + T(i):sve(i), :);
                Y(nes(i) + 2 * T(i):nee(i), :) = Y(nes(i) + 2 * T(i):nee(i), :) ...
                    + diag(silJ{i}(4 * T(i) + 2 * S(i) + 1:5 * T(i) + 2 * S(i), 1)) ...
                    * X(svs(i):sve(i) - T(i), :);
                Y(nes(i) + 2 * T(i):nee(i), :) = Y(nes(i) + 2 * T(i):nee(i), :)...
                    + diag(silJ{i}(4 * T(i) + 2 * S(i) + 1:5 * T(i) + 2 * S(i), 2)) ...
                    * X(svs(i) + T(i):sve(i), :);
                
                Y(es(n + 1) + 1:es(n + 1) + 3 * P, :) = ...
                    Y(es(n + 1) + 1:es(n + 1) + 3 * P, :) ...
                    + (parameters.xi_0 / n) * tplatesqrt_3 * ...
                    vars(vs(n + 1) + 1:vs(n + 1) + 3 * P)' * X(sv(i), :);
                for m = 1:cM
                    ms = es(n + 1) + 1 +       m * 3 * P;
                    me = es(n + 1) +     (m + 1) * 3 * P;
                    Y(ms:me, :) = Y(ms:me, :) + tplatesqrt_3 * ...
                        (parameters.xi_def / n) * vars(vs(n + 1) + 1 + m * 3 * P:...
                        vs(n + 1) + (m + 1) * 3 * P)' ...
                        * X(sv(i), :);
                end
                
                Y(mes(i):mee(i), :) = 2 * beta * X(mvs(i):mve(i), :);
            end
            
            Y(es(n + 1) + 1:es(n + 1) + 3 * P, :) = ...
                Y(es(n + 1) + 1:es(n + 1) + 3 * P, :) + ...
                parameters.xi_0 * tplatesqrt_3 * meansc * ...
                X(vs(n + 1) + 1:vs(n + 1) + 3 * P, :);
            for m = 1:cM
                ms = es(n + 1) + 1 +       m * 3 * P;
                me = es(n + 1) +     (m + 1) * 3 * P;
                
                Y(ms:me, :) = Y(ms:me, :) + parameters.xi_def * tplatesqrt_3 * ...
                    meansc * X(vs(n + 1) + 1 + m * 3 * P:...
                    vs(n + 1) + (m + 1) * 3 * P, :);
            end
            
            if both
                flag = -1;
            end
        else
            Y = X;
        end
        
        if flag < 0
            W = zeros(length(vars), size(Y, 2));
            
            for i = 1:n
                esi = es(i) + 1;
                esn = es(i + 1) - cM;
                W(rvs(i):mve(i), :) = varsJ{i}' * Y(esi:esn, :);
                
                ms  = vs(n + 1) + 1;
                me  = vs(n + 1) + 3 * P;
                W(ms:me, :) = W(ms:me, :) + meshJ{i}' * Y(esi:esn, :);
                
                for m = 1:cM
                    ms  = vs(n + 1) + 1 + m * 3 * P;
                    me  = vs(n + 1) + (m + 1) * 3 * P;
                    
                    W(ms:me, :) = W(ms:me, :) + ...
                        vars(mvs(i) + m - 1) * meshJ{i}' * Y(esi:esn, :);
                end
                
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    + diag(silJ{i}(1:T(i), 1)) * Y(ses(i):see(i) - T(i), :);
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + diag(silJ{i}(1:T(i), 2)) * Y(ses(i):see(i) - T(i), :);
                
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    + diag(silJ{i}(T(i) + 1:2 * T(i), 1)) ...
                    * Y(ses(i) + T(i):see(i), :);
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + diag(silJ{i}(T(i) + 1:2 * T(i), 2)) ...
                    * Y(ses(i) + T(i):see(i), :);
                
                D = diag(silJ{i}(2 * T(i) + 1:2 * T(i) + S(i), 1)) + ...
                    diag(silJ{i}(2 * T(i) + S(i) + 1:...
                    2 * T(i) + 2 * S(i) - 1, 1), -1);
                D(1, S(i)) = silJ{i}(2 * T(i) + 2 * S(i), 1);
                D(fixed_sil_pts{i}, :) = [];
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    +  D * Y(ces(i):cee(i), :);
                
                D = diag(silJ{i}(2 * T(i) + 1:2 * T(i) + S(i), 2)) + ...
                    diag(silJ{i}(2 * T(i) + S(i) + 1:...
                    2 * T(i) + 2 * S(i) - 1, 2), -1);
                D(1, S(i)) = silJ{i}(2 * T(i) + 2 * S(i), 2);
                D(fixed_sil_pts{i}, :) = [];
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + D * Y(ces(i):cee(i), :);
                
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    + diag(silJ{i}(2 * T(i) + 2 * S(i) + 1:3 * T(i) + 2 * S(i), 1)) ...
                    * Y(nes(i):nee(i) - 2 * T(i), :);
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + diag(silJ{i}(2 * T(i) + 2 * S(i) + 1:3 * T(i) + 2 * S(i), 2)) ...
                    * Y(nes(i):nee(i) - 2 * T(i), :);
                
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    + diag(silJ{i}(3 * T(i) + 2 * S(i) + 1:4 * T(i) + 2 * S(i), 1)) ...
                    * Y(nes(i) + T(i):nee(i) - T(i), :);
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + diag(silJ{i}(3 * T(i) + 2 * S(i) + 1:4 * T(i) + 2 * S(i), 2)) ...
                    * Y(nes(i) + T(i):nee(i) - T(i), :);
                
                W(svs(i):sve(i) - T(i), :) = W(svs(i):sve(i) - T(i), :) ...
                    + diag(silJ{i}(4 * T(i) + 2 * S(i) + 1:5 * T(i) + 2 * S(i), 1)) ...
                    * Y(nes(i) + 2 * T(i):nee(i), :);
                W(svs(i) + T(i):sve(i), :) = W(svs(i) + T(i):sve(i), :) ...
                    + diag(silJ{i}(4 * T(i) + 2 * S(i) + 1:5 * T(i) + 2 * S(i), 2)) ...
                    * Y(nes(i) + 2 * T(i):nee(i), :);
                
                W(sv(i), :) = W(sv(i), :) + (parameters.xi_0 / n) * ...
                    vars(vs(n + 1) + 1:vs(n + 1) + 3 * P) * ...
                    tplatesqrt_3' * Y(es(n + 1) + 1:es(n + 1) + 3 * P, :);
                
                for m = 1:cM
                    ms = es(n + 1) + 1 +       m * 3 * P;
                    me = es(n + 1) +     (m + 1) * 3 * P;
                    
                    W(sv(i), :) = W(sv(i), :) + (parameters.xi_def / n) * ...
                        vars(vs(n + 1) + 1 + m * 3 * P:...
                        vs(n + 1) + (m + 1) * 3 * P) * ...
                        tplatesqrt_3' * Y(ms:me, :);
                end
                
                W(mvs(i):mve(i), :) = ...
                    W(mvs(i):mve(i), :) + 2 * beta * Y(mes(i):mee(i), :);
            end
            
            W(vs(n + 1) + 1:vs(n + 1) + 3 * P, :) = ...
                W(vs(n + 1) + 1:vs(n + 1) + 3 * P, :) + ...
                parameters.xi_0 * meansc * tplatesqrt_3 * ...
                Y(es(n + 1) + 1:es(n + 1) + 3 * P, :);
            for m = 1:cM
                ms  = es(n + 1) + 1 +      m  * 3 * P;
                me  = es(n + 1) +     (m + 1) * 3 * P;
                sms = vs(n + 1) + 1 +      m  * 3 * P;
                sme = vs(n + 1) +     (m + 1) * 3 * P;
                
                W(sms:sme, :) = W(sms:sme, :) + ...
                    parameters.xi_def * tplatesqrt_3 * meansc * Y(ms:me, :);
            end
        else
            W = Y;
        end
        
        assert(~any(any(isnan(W))));
    end

% Note: To avoid a lot of 2^-0.5 factors, this function calculates 2E,
% where E is the energy described in the paper.
    function [energy Jdata] = calc_energy(vars)
        energy = zeros(es(n + 1) + 3 * P * (cM + 1), 1);
        modes  = reshape(vars(vs(n + 1) + 1:end), 3 * P, cM + 1);
        meansc = 0;
        
        for i = 1:n
            sil_barycentric{i}  = last_sil_barycen{i};
            sil_triangles{i}    = last_sil_tris{i};
            update_sil_surface(i, vars(svs(i):sve(i)));
            
            % Calculate derivatives that give the relation between the
            % parameter space for sil_barycentric (with transitions around
            % extraordinary points), (u,v), and the uniform parameter space
            % associated with each triangle, (a,b).
            dadu = zeros(S(i), 1);
            dbdu = zeros(S(i), 1);
            dadv = zeros(S(i), 1);
            dbdv = zeros(S(i), 1);
            for j = 1:S(i)
                if fixed_sil_pts{i}(j)
                    continue
                end
                
                cpt     = sil_barycentric{i}(j, 1) * u_basis + ...
                    sil_barycentric{i}(j, 2) * v_basis;
                e       = mesh.edges(sil_triangles{i}(j));
                val     = e.next.vert.valency;
                reg     = abs(cpt) > vertex_radius || val == 6;
                
                if reg
                    dadu(j) = 1;
                    dbdu(j) = 0;
                    dadv(j) = -1 / sqrt(3);
                    dbdv(j) = 2 / sqrt(3);
                elseif cpt == 0
                    dadu(j) = 1;
                    dbdu(j) = 0;
                    if val == 3
                        dbdv(j) = cos(atan(sqrt(3) / 2));
                        dadv(j) = 0.5 * dbdv(j);
                    elseif val == 4
                        dadv(j) = 0;
                        dbdv(j) = 1;
                    else
                        dv      = [ -1/sin(pi/val) - 1/cos(2*pi/val) ...
                            -1/sin(pi/val) ];
                        dv      = dv ./ abs(-dv(1) * v_basis + dv(2));
                        dadv(j) = dv(2);
                        dbdv(j) = -dv(1);
                    end
                else
                    cpt     = cpt ^ (6 / val);
                    re      = real(cpt);
                    im      = imag(cpt);
                    lsq     = re^2 + im^2;
                    theta   = val * angle(cpt) / 6;
                    co      = cos(theta);
                    sn      = sin(theta);
                    r3      = sqrt(3);
                    dadu(j) = (val * lsq ^ (val / 12 - 1)) * ...
                        (     re * co / 6  +      im * sn / 6  + ...
                        r3 * im * co / 18 - r3 * re * sn / 18);
                    dadv(j) = (-val * lsq ^ (val / 12 - 1)) * ...
                        (    -im * co / 6  +      re * sn / 6  + ...
                        r3 * re * co / 18 + r3 * im * sn / 18);
                    dbdu(j) = -(r3 * val * lsq ^ (val / 12 - 1) * ...
                        (im * co - re * sn)) / 9;
                    dbdv(j) =  (r3 * val * lsq ^ (val / 12 - 1) * ...
                        (re * co + im * sn)) / 9;
                end
            end
            
            v   = vars(rvs(i):rve(i));
            if nargout > 1
                silJ{i}  = zeros(2 * S(i) + 5 * T(i),            2     );
                varsJ{i} = zeros(    S(i) + 5 * T(i) + 2 * K(i), 7 + cM);
                meshJ{i} = zeros(    S(i) + 5 * T(i) + 2 * K(i), 3 * P );
                
                % To calculate rotation derivatives
                vlen         = norm(v, 2);
                hsinc        = sinc_unnorm(vlen / 2) / 2;
                dqdv         = zeros(4, 3);
                dqdv(4, :)   = -0.5 * v * hsinc;
                if vlen < eps^(1/4)
                    % Use a Taylor expansion approximation to avoid
                    % numerical instability
                    mult     = ((vlen ^ 2 / 40) - 1) / 24;
                else
                    mult     = (cos(vlen / 2) / 2 - hsinc) / (vlen ^ 2);
                end
                dqdv(1:3, :) = repmat(v, 1, 3) .* repmat(v', 3, 1) .* mult;
                dqdv(1, 1)   = dqdv(1, 1) + hsinc;
                dqdv(2, 2)   = dqdv(2, 2) + hsinc;
                dqdv(3, 3)   = dqdv(3, 3) + hsinc;
                
                q            = [ hsinc * v ; cos(vlen / 2) ];
                
                dRdq         = 2 * ...
                    [         0  -2 * q(2)  -2 * q(3)      0 ;
                    q(2)       q(1)       q(4)   q(3) ;
                    q(3)      -q(4)       q(1)  -q(2) ;
                    q(2)       q(1)      -q(4)  -q(3) ;
                    -2 * q(1)          0  -2 * q(3)      0 ;
                    q(4)       q(3)       q(2)   q(1) ;
                    q(3)       q(4)       q(1)   q(2) ;
                    -q(4)       q(3)       q(2)  -q(1) ;
                    -2 * q(1)  -2 * q(2)          0      0 ];
                
                dRdv         = dRdq * dqdv;
                rs           = project.images(i).transform;
                dMdv_noscale = [ rs(1:3, 1:3) * dRdv(1:3, :) ;
                    rs(1:3, 1:3) * dRdv(4:6, :) ;
                    rs(1:3, 1:3) * dRdv(7:9, :) ];
                dMdv         = dMdv_noscale * vars(sv(i));
            end
            
            % Silhouette points
            rot         = [ expm([    0 -v(3)  v(2) ;
                v(3)     0 -v(1) ;
                -v(2)  v(1)     0 ]) zeros(3, 1) ;
                zeros(1, 3) 1 ];
            rotscale    = rot .* vars(sv(i));
            rotscale(4, 4) = 1;
            translate   = [ 1 0 0 vars(tvs(i)    ) ;
                0 1 0 vars(tvs(i) + 1) ;
                0 0 1 vars(tvs(i) + 2) ;
                0 0 0 1                ];
            DT          = translate * ...
                project.images(i).transform * rotscale;
            noscale     = translate * project.images(i).transform * rot;
            pointvars   = reshape(modes * [1 ; vars(mvs(i):mve(i))], P, 3);
            transformed = [ pointvars ones(P, 1) ] * DT';
            transformed = transformed(:, 1:3);
            rotscaled   = [ pointvars ones(P, 1) ] * rotscale' * ...
                project.images(i).transform';
            rotscaled   = rotscaled(:, 1:3);
            rotated     = [ pointvars ones(P, 1) ] * rot' * ...
                project.images(i).transform';
            rotated     = rotated(:, 1:3);
            jT          = 0;
            for j = 0:S(i) - 1
                % Continuity
                if j == S(i) - 1, nextj = 1; else nextj = j + 2; end
                [dist ddda1 dddb1 ddda2 dddb2] = mesh.distbetween(...
                    mesh.edges(sil_triangles{i}(j + 1)), ...
                    sil_barycentric{i}(j + 1, 1), ...
                    sil_barycentric{i}(j + 1, 2), ...
                    mesh.edges(sil_triangles{i}(nextj)), ...
                    sil_barycentric{i}(nextj, 1), ...
                    sil_barycentric{i}(nextj, 2));
                energy(ces(i) + j) = sqrt(2 * gamma) * dist;
                
                if nargout > 1
                    dddu1 = ddda1 * dadu(j + 1) + dddb1 * dbdu(j + 1);
                    dddv1 = ddda1 * dadv(j + 1) + dddb1 * dbdv(j + 1);
                    dddu2 = ddda2 * dadu(nextj) + dddb2 * dbdu(nextj);
                    dddv2 = ddda2 * dadv(nextj) + dddb2 * dbdv(nextj);
                    rt2g  = sqrt(2 * gamma);
                    
                    silJ{i}(2 * T(i)        + j + 1, 1) = rt2g * dddu1;
                    silJ{i}(2 * T(i)        + j + 1, 2) = rt2g * dddv1;
                    silJ{i}(2 * T(i) + S(i) + j + 1, 1) = rt2g * dddu2;
                    silJ{i}(2 * T(i) + S(i) + j + 1, 2) = rt2g * dddv2;
                end
                
                if fixed_sil_pts{i}(j + 1)
                    continue;
                end
                
                [lt deriv sec] = mesh.limitevaluation(...
                    mesh.edges(sil_triangles{i}(j + 1)), ...
                    sil_barycentric{i}(j + 1, 1), ...
                    sil_barycentric{i}(j + 1, 2));
                % X
                dx = lt * transformed(:, 1) - sil_pts{i}(j + 1, 1);
                
                % Y
                dy = lt * transformed(:, 2) - sil_pts{i}(j + 1, 2);
                
                energy(ses(i)        + jT) = dx / parameters.sigma_sil;
                energy(ses(i) + T(i) + jT) = dy / parameters.sigma_sil;
                
                if nargout > 1
                    % X
                    xd = lt' * DT(1, 1:3);
                    xd = xd(:);
                    meshJ{i}(       jT + 1, :) = xd / parameters.sigma_sil;
                    
                    % Y
                    yd = lt' * DT(2, 1:3);
                    yd = yd(:);
                    meshJ{i}(T(i) + jT + 1, :) = yd / parameters.sigma_sil;
                end
                
                % U and V
                dxda    = [  0 1 0 ] * deriv;
                dxdb    = [ -1 0 0 ] * deriv;
                dxdu    = dxda * dadu(j + 1) + dxdb * dbdu(j + 1);
                dxdv    = dxda * dadv(j + 1) + dxdb * dbdv(j + 1);
                
                if nargout > 1
                    dxduval                      = dxdu * rotscaled;
                    dxdvval                      = dxdv * rotscaled;
                    silJ{i}(jT + 1, 1)          = dxduval(1);
                    silJ{i}(jT + 1, 2)          = dxdvval(1);
                    varsJ{i}(jT + 1, 1:3)        = lt * pointvars * ...
                        dMdv([1 4 7], :);
                    varsJ{i}(jT + 1, 4)          = 1;
                    varsJ{i}(jT + 1, 7)          = lt * pointvars * ...
                        noscale(1, 1:3)';
                    
                    silJ{i}(T(i) + jT + 1, 1)   = dxduval(2);
                    silJ{i}(T(i) + jT + 1, 2)   = dxdvval(2);
                    varsJ{i}(T(i) + jT + 1, 1:3) = lt * pointvars * ...
                        dMdv([2 5 8], :);
                    varsJ{i}(T(i) + jT + 1, 5)   = 1;
                    varsJ{i}(T(i) + jT + 1, 7)   = lt * pointvars * ...
                        noscale(2, 1:3)';
                    
                    for m = 1:cM
                        varsJ{i}(       jT + 1, 7 + m) = ...
                            meshJ{i}(       jT + 1, :) * modes(:, m + 1);
                        varsJ{i}(T(i) + jT + 1, 7 + m) = ...
                            meshJ{i}(T(i) + jT + 1, :) * modes(:, m + 1);
                    end
                end
                
                % Normals
                dxduval    = dxdu * rotated;
                dxdvval    = dxdv * rotated;
                unormal    = cross(dxdvval, dxduval);
                normlen    = sqrt(sum(unormal .^ 2));
                normal     = unormal / normlen;
                s_norm_rep = 1 / parameters.sigma_norm;
                
                energy(nes(i) +            jT) = s_norm_rep * (normal(1) - ...
                    sil_normals{i}(j + 1, 1));
                energy(nes(i) +     T(i) + jT) = s_norm_rep * (normal(2) - ...
                    sil_normals{i}(j + 1, 2));
                energy(nes(i) + 2 * T(i) + jT) = s_norm_rep * normal(3);
                
                if nargout > 1
                    cpt     = sil_barycentric{i}(j + 1, 1) * u_basis + ...
                        sil_barycentric{i}(j + 1, 2) * v_basis;
                    e       = mesh.edges(sil_triangles{i}(j + 1));
                    val     = e.next.vert.valency;
                    reg     = abs(cpt) > vertex_radius || val == 6;
                    if reg
                        d2xdu2  = sec(1, :);
                        d2xdudv = ([ -1 2 0 ] / sqrt(3)) * sec;
                        d2xdv2  = ([ 1 -4 4 ] / 3) * sec;
                    elseif cpt == 0
                        if val == 3
                            dv   = [ -1 0.5 0 ] * cos(atan(sqrt(3) / 2));
                        elseif val == 4
                            dv   = [ -1 0 0 ];
                        else
                            dv   = [ -1/sin(pi / val) - 1/cos(2 * pi / val) ...
                                -1/sin(pi / val) 0 ];
                            dv   = dv ./ abs(-dv(1) * v_basis + dv(2) + ...
                                dv(3) * (v_basis - u_basis));
                        end
                        
                        d2xdu2  = sec(1, :);
                        d2xdudv = [ dv(2) -dv(1) 0 ] * sec;
                        d2xdv2  = [ dv(2) ^ 2  -2 * dv(1) * dv(2) ...
                            dv(1) ^ 2 ] * sec;
                    else
                        cpt   = cpt ^ (6 / val);
                        re    = real(cpt);
                        im    = imag(cpt);
                        lsq   = re^2 + im^2;
                        theta = val * angle(cpt) / 6;
                        co    = cos(theta);
                        sn    = sin(theta);
                        r3    = sqrt(3);
                        rere  = re ^ 2;
                        imim  = im ^ 2;
                        reim  = re * im;
                        
                        d2adu2  = val * lsq ^ (val / 12 - 2) * ...
                            (val - 6) * (     rere * co / 36  - ...
                            imim * co / 36  + ...
                            reim * sn / 18  - ...
                            r3 * rere * sn / 108 + ...
                            r3 * imim * sn / 108 + ...
                            r3 * reim * co / 54);
                        d2adudv = -val * lsq ^ (val / 12 - 2) * ...
                            (val - 6) * (     rere * sn / 36  - ...
                            imim * sn / 36  - ...
                            reim * co / 18  + ...
                            r3 * rere * co / 108 - ...
                            r3 * imim * co / 108 + ...
                            r3 * reim * sn / 54);
                        d2adv2  = -d2adu2;
                        
                        d2bdu2  = -(r3 * val * lsq ^ (val / 12 - 2) * ...
                            (val - 6) * (-rere * sn + ...
                            imim * sn + 2 * reim * co)) / 54;
                        d2bdudv =  (r3 * val * lsq ^ (val / 12 - 2) * ...
                            (val - 6) * ( rere * co - ...
                            imim * co + 2 * reim * sn)) / 54;
                        d2bdv2  = -d2bdu2;
                        
                        d2xdu2  = [ dadu(j + 1)^2 ...
                            2 * dadu(j + 1) * dbdu(j + 1) ...
                            dbdu(j + 1)^2 ] * sec + ...
                            d2adu2 * dxda + d2bdu2 * dxdb;
                        d2xdudv = [ dadu(j + 1) * dadv(j + 1) ...
                            dadu(j + 1) * dbdv(j + 1) + ...
                            dbdu(j + 1) * dadv(j + 1) ...
                            dbdu(j + 1) * dbdv(j + 1) ] * sec + ...
                            d2adudv * dxda + d2bdudv * dxdb;
                        d2xdv2  = [ dadv(j + 1)^2 ...
                            2 * dadv(j + 1) * dbdv(j + 1) ...
                            dbdv(j + 1)^2 ] * sec + ...
                            d2adv2 * dxda + d2bdv2 * dxdb;
                    end
                    
                    d2xdu2  = d2xdu2  * rotated;
                    d2xdudv = d2xdudv * rotated;
                    d2xdv2  = d2xdv2  * rotated;
                    
                    % Changes to the mesh
                    dundx   = [ zeros(1, P)                           ;
                        dxdu * dxdvval(3) - dxdv * dxduval(3) ;
                        dxdv * dxduval(2) - dxdu * dxdvval(2) ];
                    dundy   = [ dxdv * dxduval(3) - dxdu * dxdvval(3) ;
                        zeros(1, P)                           ;
                        dxdu * dxdvval(1) - dxdv * dxduval(1) ];
                    dundz   = [ dxdu * dxdvval(2) - dxdv * dxduval(2) ;
                        dxdv * dxduval(1) - dxdu * dxdvval(1) ;
                        zeros(1, P)                           ];
                    
                    dlndx   = normal * dundx;
                    dlndy   = normal * dundy;
                    dlndz   = normal * dundz;
                    
                    dndx    = (dundx - normal' * dlndx) / normlen;
                    dndy    = (dundy - normal' * dlndy) / normlen;
                    dndz    = (dundz - normal' * dlndz) / normlen;
                    
                    n1ej    = S(i) + 2 * T(i) + jT + 1;
                    n2ej    = S(i) + 3 * T(i) + jT + 1;
                    n3ej    = S(i) + 4 * T(i) + jT + 1;
                    
                    meshJ{i}(n1ej, 1:P)             = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 1) + ...
                        dndy(1,:) * noscale(2, 1) + ...
                        dndz(1,:) * noscale(3, 1));
                    meshJ{i}(n1ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 2) + ...
                        dndy(1,:) * noscale(2, 2) + ...
                        dndz(1,:) * noscale(3, 2));
                    meshJ{i}(n1ej, 2 * P + 1:3 * P) = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 3) + ...
                        dndy(1,:) * noscale(2, 3) + ...
                        dndz(1,:) * noscale(3, 3));
                    
                    meshJ{i}(n2ej, 1:P)             = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 1) + ...
                        dndy(2,:) * noscale(2, 1) + ...
                        dndz(2,:) * noscale(3, 1));
                    meshJ{i}(n2ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 2) + ...
                        dndy(2,:) * noscale(2, 2) + ...
                        dndz(2,:) * noscale(3, 2));
                    meshJ{i}(n2ej, 2 * P + 1:3 * P) = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 3) + ...
                        dndy(2,:) * noscale(2, 3) + ...
                        dndz(2,:) * noscale(3, 3));
                    
                    meshJ{i}(n3ej, 1:P)             = s_norm_rep * ...
                        (dndx(3,:) * noscale(1, 1) + ...
                        dndy(3,:) * noscale(2, 1) + ...
                        dndz(3,:) * noscale(3, 1));
                    
                    meshJ{i}(n3ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(3,:) * noscale(1, 2) + ...
                        dndy(3,:) * noscale(2, 2) + ...
                        dndz(3,:) * noscale(3, 2));
                    
                    meshJ{i}(n3ej, 2 * P + 1:3 * P) = s_norm_rep * ...
                        (dndx(3,:) * noscale(1, 3) + ...
                        dndy(3,:) * noscale(2, 3) + ...
                        dndz(3,:) * noscale(3, 3));
                    
                    
                    % Changes to the silhouette points
                    dundu = (cross(dxdvval, d2xdu2 ) + ...
                        cross(d2xdudv, dxduval))';
                    dundv = (cross(dxdvval, d2xdudv) + ...
                        cross(d2xdv2 , dxduval))';
                    
                    dlndu = normal * dundu;
                    dlndv = normal * dundv;
                    
                    dndu  = (dundu - normal' * dlndu) / normlen;
                    dndv  = (dundv - normal' * dlndv) / normlen;
                    
                    silJ{i}(n1ej + S(i), 1) = s_norm_rep * dndu(1);
                    silJ{i}(n1ej + S(i), 2) = s_norm_rep * dndv(1);
                    silJ{i}(n2ej + S(i), 1) = s_norm_rep * dndu(2);
                    silJ{i}(n2ej + S(i), 2) = s_norm_rep * dndv(2);
                    silJ{i}(n3ej + S(i), 1) = s_norm_rep * dndu(3);
                    silJ{i}(n3ej + S(i), 2) = s_norm_rep * dndv(3);
                    
                    % Changes to the matrix transform from rotation
                    dndvar = ...
                        (dndx * pointvars * dMdv_noscale([1 4 7], :) + ...
                        dndy * pointvars * dMdv_noscale([2 5 8], :) + ...
                        dndz * pointvars * dMdv_noscale([3 6 9], :));
                    
                    varsJ{i}(n1ej, 1:3) = s_norm_rep * dndvar(1, :);
                    varsJ{i}(n2ej, 1:3) = s_norm_rep * dndvar(2, :);
                    varsJ{i}(n3ej, 1:3) = s_norm_rep * dndvar(3, :);
                    
                    for m = 1:cM
                        varsJ{i}(n1ej, 7 + m) = ...
                            meshJ{i}(n1ej, :) * modes(:, m + 1);
                        varsJ{i}(n2ej, 7 + m) = ...
                            meshJ{i}(n2ej, :) * modes(:, m + 1);
                        varsJ{i}(n3ej, 7 + m) = ...
                            meshJ{i}(n3ej, :) * modes(:, m + 1);
                    end
                end
                
                jT = jT + 1;
            end
            
            % Constraint vertices
            energy(cxes(i):cyee(i)) = repmat(conwt{i}, 2, 1) .* ...
                (constAeq{i} * reshape(transformed(:, 1:2), 2 * P, 1) - ...
                constbeq{i});
            if nargout > 1
                C = repmat(conwt{i}, 1, P) .* constAeq{i}(1:K(i), 1:P);
                
                % Derivatives for equations measuring difference in X
                cxs = S(i) + 5 * T(i) + 1;
                cxe = S(i) + 5 * T(i) + K(i);
                % - Changes to mesh
                meshJ{i}(cxs:cxe,         1:P)     = DT(1, 1) * C;
                meshJ{i}(cxs:cxe,     P + 1:2 * P) = DT(1, 2) * C;
                meshJ{i}(cxs:cxe, 2 * P + 1:3 * P) = DT(1, 3) * C;
                % - Changes to matrix transform
                %   - From translation
                varsJ{i}(cxs:cxe, 4)               = conwt{i};
                %   - From rotation
                varsJ{i}(cxs:cxe, 1:3)             = ...
                    C * pointvars * dMdv([1 4 7], :);
                %   - From scaling
                varsJ{i}(cxs:cxe, 7)               = ...
                    C * pointvars * noscale(1, 1:3)';
                
                % Derivatives for equations measuring difference in Y
                cys = S(i) + 5 * T(i) +     K(i) + 1;
                cye = S(i) + 5 * T(i) + 2 * K(i);
                % - Changes to mesh
                meshJ{i}(cys:cye,         1:P)     = DT(2, 1) * C;
                meshJ{i}(cys:cye,     P + 1:2 * P) = DT(2, 2) * C;
                meshJ{i}(cys:cye, 2 * P + 1:3 * P) = DT(2, 3) * C;
                % - Changes to matrix transform
                %   - From translation
                varsJ{i}(cys:cye, 5)               = conwt{i};
                %   - From rotation
                varsJ{i}(cys:cye, 1:3)             = ...
                    C * pointvars * dMdv([2 5 8], :);
                %   - From scaling
                varsJ{i}(cys:cye, 7)               = ...
                    C * pointvars * noscale(2, 1:3)';
                
                for m = 1:cM
                    varsJ{i}(cxs:cye, 7 + m) = ...
                        meshJ{i}(cxs:cye, :) * modes(:, m + 1);
                end
            end
            
            % Shape coefficient regularisation
            energy(mes(i):mee(i)) = 2 * beta * vars(mvs(i):mve(i));
            meansc                = meansc + vars(sv(i));
        end
        meansc = meansc / n;
        
        % Shape mode regularisation
        energy(es(n + 1) + 1:es(n + 1) + 3 * P) = ...
            parameters.xi_0 * meansc * tplatesqrt_3 * modes(:, 1);
        for m = 1:cM
            ms = es(n + 1) + 1 +       m * 3 * P;
            me = es(n + 1) +     (m + 1) * 3 * P;
            energy(ms:me) = parameters.xi_def * meansc * tplatesqrt_3 * modes(:, m + 1);
        end
        
        if nargout > 1
            assert(max(S) > 2 * max(K));
            Jdata = zeros(2 * max(S) + 5 * max(T) + 1, ...
                n * (9 + cM + 3 * P));
            Jdata(end, 1:length(vars)) = vars';
            
            for i = 1:n
                bs = (i - 1) * (9 + cM + 3 * P) + 1;
                be =       i * (9 + cM + 3 * P);
                Jdata(1:2 * S(i) + 5 * T(i), bs:bs + 1) = silJ{i};
                Jdata(1:    S(i) + 5 * T(i) + 2 * K(i), ...
                    bs + 2:bs + 8 + cM) = varsJ{i};
                Jdata(1:    S(i) + 5 * T(i) + 2 * K(i), ...
                    bs + 9 + cM:be)     = meshJ{i};
            end
            
            if ~jacobmult_on
                JdataT = jacobmult(Jdata, eye(length(energy)), -1);
                Jdata  = jacobmult(Jdata, eye(length(vars)),    1);
                assert(~any(any(abs(JdataT' - Jdata) > 1e-10)));
            end
        end
    end

    function update_sil_surface(im_no, new_uvs)
        if all(last_sil_uvs{im_no} == new_uvs)
            return
        end
        
        jT = 1;
        
        for j = 1:S(im_no)
            if fixed_sil_pts{im_no}(j)
                continue;
            else
                diff    =             new_uvs([jT T(im_no) + jT]) - ...
                    last_sil_uvs{im_no}([jT T(im_no) + jT]);
                diff    = complex(diff(1), diff(2));
                jT = jT + 1;
                
                if diff == 0
                    continue;
                end
            end
            
            e           = mesh.edges(sil_triangles{im_no}(j));
            last_sil_pt = sil_barycentric{im_no}(j, 1) * u_basis + ...
                sil_barycentric{im_no}(j, 2) * v_basis;
            done        = false;
            val         = e.next.vert.valency;
            outside_crc = abs(last_sil_pt) > vertex_radius;
            
            while ~done
                % alpha is the proportion of 'diff' that we're going to
                % move in this step
                alpha     = 1;
                
                % moverot is the number of rotations around e.next.vert
                % that are required at the end of this step
                moverot   = 0;
                
                % Whether this transition occurs entirely outside a
                % 'vertex circle'
                out_trans = false;
                
                if outside_crc
                    % Rotate around the triangle so that 'diff' is towards
                    % the origin
                    diffangle = angle(diff * exp(complex(0, -pi / 6)));
                    if diffangle < 2 * pi / 3 && diffangle > -2 * pi / 3
                        u = sil_barycentric{im_no}(j, 1);
                        v = sil_barycentric{im_no}(j, 2);
                        w = 1 - (u + v);
                        if diffangle > 0
                            e           = e.next;
                            last_sil_pt = w * u_basis + u * v_basis;
                            diff        = diff * ...
                                exp(complex(0, 2 * pi / 3));
                        else
                            e           = e.next.next;
                            last_sil_pt = v * u_basis + w * v_basis;
                            diff        = diff * ...
                                exp(complex(0, -2 * pi / 3));
                        end
                        val = e.next.vert.valency;
                    end
                    
                    % Where we would like to move to, if we could
                    moveto = last_sil_pt + diff;
                    
                    if abs(moveto) < vertex_radius || ...
                            angle(moveto) > pi / 3 || angle(moveto) < 0
                        % 'moveto' is too far away: we'll need to take an
                        % interim step and keep moving.
                        
                        % Now we have to figure out whether the first
                        % intersection is with the 'vertex circle', or with
                        % the edges of the current triangle. I'll use an
                        % approximation here based on the angle of diff.
                        % Which might mean that we miss very thin snicks of
                        % the circle surrounding a vertex. But those
                        % shouldn't make much difference anyway, so this is
                        % good for efficiency.
                        top_ang = angle(circ_top - last_sil_pt);
                        if top_ang < 0
                            top_ang = top_ang + 2 * pi;
                        end
                        
                        bot_ang = angle(circ_bot - last_sil_pt);
                        if bot_ang < 0
                            bot_ang = bot_ang + 2 * pi;
                        end
                        
                        diff_ang = angle(diff);
                        if diff_ang < 0
                            diff_ang = diff_ang + 2 * pi;
                        end
                        
                        if diff_ang < top_ang
                            % Intersection with the v_basis edge
                            alpha     = ...
                                (real(v_basis) * imag(last_sil_pt) - ...
                                real(last_sil_pt) * imag(v_basis)) / ...
                                (real(diff) * imag(v_basis) - ...
                                real(v_basis) * imag(diff));
                            moverot   = 1;
                            out_trans = true;
                        elseif diff_ang > bot_ang
                            % Intersection with the u_basis edge
                            alpha     = -imag(last_sil_pt) / imag(diff);
                            moverot   = -1;
                            out_trans = true;
                        else
                            % Moving into the vertex circle
                            quadratic   = [ real(diff)^2 + imag(diff)^2
                                2 * real(last_sil_pt * ...
                                conj(diff))
                                real(last_sil_pt) ^ 2 + ...
                                imag(last_sil_pt) ^ 2 - ...
                                vertex_radius ^ 2 ];
                            alpha       = real(min(roots(quadratic)));
                        end
                    else
                        % moveto is the new position and we're done
                        done        = true;
                    end
                else
                    % We're inside the vertex circle
                    last_sil_pt = last_sil_pt   ^ (6 / val);
                    circle_r    = vertex_radius ^ (6 / val);
                    
                    % Where we would like to move to, if we could
                    moveto      = last_sil_pt + diff;
                    
                    if abs(moveto) <= circle_r
                        % moveto is the new position and we're done
                        done   = true;
                    else
                        % 'moveto' is too far away: we'll need to take an
                        % interim step and keep moving.
                        quadratic = [ real(diff)^2 + imag(diff)^2
                            2 * real(last_sil_pt * conj(diff))
                            real(last_sil_pt) ^ 2 + ...
                            imag(last_sil_pt) ^ 2 - ...
                            circle_r ^ 2 ];
                        alpha       = real(max(roots(quadratic)));
                    end
                end
                
                % Where we're actually moving to
                moveto = last_sil_pt + alpha * diff;
                diff   = (1 - alpha) * diff;
                
                % Now we rotate around e.next.vert as required
                if out_trans
                    rotvector = exp(complex(0, -pi * moverot / 3));
                elseif ~outside_crc
                    moverot   = floor((val / 2) * (angle(moveto) / pi));
                    rotvector = exp(complex(0, -2 * pi * moverot / val));
                else
                    rotvector = 1;
                end
                diff   = diff   * rotvector;
                moveto = moveto * rotvector;
                
                if ~out_trans
                    if outside_crc
                        % If we're entering a vertex circle, we need to
                        % correct the angle of 'diff' so that it makes
                        % sense in the z ^ (6 / valency) space.
                        last_sil_pt = last_sil_pt * rotvector;
                        newdiffang  = angle(moveto      ^ (6 / val) ...
                            - last_sil_pt ^ (6 / val));
                        diff = abs(diff) * exp(complex(0, newdiffang));
                    else
                        % If we're exiting a vertex circle, we need to map
                        % the extraordinary sector back into a 6-valent
                        % sector
                        mvtang = angle(moveto);
                        moveto = moveto ^ (val / 6);
                        
                        % And rotate 'diff' accordingly
                        diff   = diff * exp(complex(0, ...
                            angle(moveto) - mvtang));
                    end
                    
                    % If it's not an 'outside transition', then we've
                    % crossed into or out of a vertex circle.
                    outside_crc = ~outside_crc;
                end
                
                % For the next iteration, if there is one.
                last_sil_pt = moveto;
                
                % Rotate around e.next.vert according to moverot
                if moverot > 0
                    e = e.next.next;
                    while moverot > 0
                        e = e.pair.next;
                        moverot = moverot - 1;
                    end
                    e = e.next;
                elseif moverot < 0
                    e = e.next;
                    while moverot < 0
                        e = e.pair.next.next;
                        moverot = moverot + 1;
                    end
                    e = e.next.next;
                end
            end
            
            % New position on the surface
            sil_barycentric{im_no}(j, :) = ...
                [ real(moveto) imag(moveto) ] / change_basis';
            w = 1 - sum(sil_barycentric{im_no}(j, :));
            [~, tripos] = max([sil_barycentric{im_no}(j, :) w]);
            if tripos == 1
                sil_barycentric{im_no}(j, 1) = ...
                    sil_barycentric{im_no}(j, 2);
                sil_barycentric{im_no}(j, 2) = w;
                e = e.next.next;
            elseif tripos == 2
                sil_barycentric{im_no}(j, 2) = ...
                    sil_barycentric{im_no}(j, 1);
                sil_barycentric{im_no}(j, 1) = w;
                e = e.next;
            end
            sil_triangles{im_no}(j) = e.index_in_mesh;
        end
    end

% --
% -- M-Lint options for this file
% --
%#ok<*DEFNU>
%#ok<*FXUP>
end
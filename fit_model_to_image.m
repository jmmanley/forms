function [shapevars,result,resultVec,init_transformed,transformed,sil_pts,init_pointvars] = fit_model_to_image(project, model, image)

% Jason Manley, Oct 2017

% User-defined connectivity (mesh)
mesh             = project.mesh;
% Thin plate energy
thinplate        = mesh.thinplate();
tplatesqrt       = real(sqrtm(thinplate / 2));
tplatesqrt_3     = blkdiag(tplatesqrt, tplatesqrt, tplatesqrt);
% Set this to false for slow, checked Jacobian computation
jacobmult_on     = true;

M = model.M;                                    % # dimensions in model
cM = M;                                         % # dimensions in calc_energy
S = 125;                                        % resolution of silhouette
T = S - sum(image.constraintsonsil);            % # free contour generator points
K = length(image.constraints3d);                % # constraint points
conwt = ones(K,1) / model.parameters.sigma_con; % ???
P = length(project.vertices);                   % # vertices in mesh model
if isfield(model.parameters,'beta')
    beta = model.parameters.beta;
else
    beta = 0.5;
end
if isfield(model.parameters,'gamma')
    gamma = model.parameters.gamma;
else
    gamma = 1/128;
end

% basis shapes
modes = model.shapemodes;

% --
% -- Constraint points
% --

constverts      = image.constraints3d;
constAeq     = zeros(2 * K, 2 * P);

for k = 1:K
    e           = mesh.verts(constverts(k)).edge.next;
    limitpoint  = mesh.limitevaluation(e, 0, 0);
    
    % X
    constAeq(    k,     1:    P) = limitpoint;
    % Y
    constAeq(K + k, P + 1:2 * P) = limitpoint;
end
constbeq     = image.constraints2d(:);



% --
% -- A large pile of index vectors to make life easier
% --

vs(1)  = 0;
vs(2)  = vs(1) + 2 * T + 7 + cM;
svs    = vs(1) + 1;
sve    = vs(1) + 2 * T;
rvs    = vs(1) + 2 * T + 1;
rve    = vs(1) + 2 * T + 3;
tvs    = vs(1) + 2 * T + 4;
tve    = vs(1) + 2 * T + 6;
sv     = vs(1) + 2 * T + 7;
mvs    = vs(1) + 2 * T + 8;
mve    = vs(1) + 2 * T + 7 + cM;

es(1) = 0;
es(2)  = es(1) + 5 * T + S + 2 * K + cM;
ses    = es(1) + 1;
see    = es(1) + 2 * T;
ces    = es(1) + 2 * T + 1;
cee    = es(1) + 2 * T + S;
nes    = es(1) + 2 * T + S + 1;
nee    = es(1) + 5 * T + S;
cxes   = es(1) + 5 * T + S + 1;
cxee   = es(1) + 5 * T + S +     K;
cyes   = es(1) + 5 * T + S +     K + 1;
cyee   = es(1) + 5 * T + S + 2 * K;
mes    = es(1) + 5 * T + S + 2 * K + 1;
mee    = es(1) + 5 * T + S + 2 * K + cM;


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
% -- Initialization for the optimizer
% --

last_sil_uvs = zeros(2 * T, 1);
problem.x0          = zeros(vs(2) + (cM + 1) * 3 * P, 1);
% problem.x0(1:vs(2)-cM) = [last_sil_uvs; zeros(6,1); 1];
% problem.x0 = [problem.x0; project.vertices(:)];
% problem.x0 = [problem.x0; zeros(3 * P, 1)];
silJ                = [];
varsJ               = [];
meshJ               = [];


% --
% -- Silhouette points
% --

fixed_sil_pts     = false(S, 1);
silhouette        = image.points;
sil_params        = sil_sample(S, silhouette);
sil_pts           = zeros(S, 2);
sil_normals       = zeros(S, 2);

for s = 1:S
    zerotangent = true;
    while zerotangent
        seg           = floor(sil_params(s));
        sil_pts(s, :) = sil_evalbezier(...
            silhouette(:, :, 1 + seg), ...
            sil_params(s) - seg);
        tan_pts       = 3 * (silhouette(2:end, :, 1 + seg) - ...
            silhouette(1:end - 1, :, 1 + seg));
        tangent       = sil_evalbezier(tan_pts, ...
            sil_params(s) - seg);
        
        zerotangent   = (norm(tangent, 2) == 0);
        % Next place to try sampling the silhouette, if it turned
        % out that the current place has zero first derivative.
        % (sil_params is not used again, so it's safe to modify
        % it).
        sil_params(s) = sil_params(s) + 1e-4;
    end
    sil_tan = tangent / norm(tangent, 2);
    
    if image.normalsLeft
        sil_normals(s, :) = [ -sil_tan(2)  sil_tan(1) ];
    else
        sil_normals(s, :) = [  sil_tan(2) -sil_tan(1) ];
    end
end

mu            = zeros(sum(image.constraintsonsil), 1);
fixed_indices = zeros(sum(image.constraintsonsil), 1);
f  = 0;
for k = 1:K
    if image.constraintsonsil(k)
        % Snap constraint points on silhouette to nearest
        % silhouette point
        closest = sil_pts - ...
            repmat(image.constraints2d(k, :), S, 1);
        [~, ind]           = min(sum(closest .^ 2, 2));
        constbeq(k)        = sil_pts(ind, 1);
        constbeq(K + k) = sil_pts(ind, 2);
        assert(fixed_sil_pts(ind) == false);
        fixed_sil_pts(ind) = true;
        f                  = f + 1;
        mu(f)              = image.constraints3d(k);
        fixed_indices(f)   = ind;
    end
end

[~, correct_mu] = sort(fixed_indices);
mu = mu(correct_mu);


init_sigma_sil   = model.parameters.sigma_sil;
init_sigma_norm  = model.parameters.sigma_norm;
init_gamma       = gamma;

fprintf(1, '\nFinding contour generator\n');
problem.x0(vs(2) + 1:end) = reshape(modes,size(problem.x0(vs(2) + 1:end)));

problem.x0(sv)   = 1;
DT               = image.transform;

pointvars        = reshape(modes * ...
    [1 ; problem.x0(mvs:mve)], ...
    P, 3);
init_pointvars = pointvars;
transformed      = [ pointvars ones(P, 1) ] * DT';
transformed      = transformed(:, 1:3);
init_transformed = transformed;
cand_norms       = cross(project.cand_derivs(:, :, 2) * ...
    transformed, ...
    project.cand_derivs(:, :, 1) * ...
    transformed);
cand_norms       = cand_norms ./ ...
    repmat(sqrt(sum(cand_norms .^ 2, 2)), 1, 3);
transformed      = transformed(:, 1:2);

if any(fixed_sil_pts)
    sil_selcands = zeros(1, S);
    nf           = sum(fixed_sil_pts);
    fixed_cands  = find(fixed_sil_pts);
    for f = 1:nf
        if f == nf, g = 1; else g = f + 1; end
        
        lin_path = sil_conspreimage(sil_pts, ...
            sil_normals, init_sigma_sil, ...
            init_gamma, init_sigma_norm, ...
            project.cand_limits * transformed, ...
            project.cand_dists, cand_norms, ...
            [mu(f) mu(f)], [mu(g) mu(g)], ...
            [fixed_cands(f) fixed_cands(g)], false);
        
        if f < nf
            sil_selcands(fixed_cands(f):fixed_cands(g)) = ...
                lin_path;
        else
            sil_selcands(fixed_cands(f):end) = ...
                lin_path(1:S - fixed_cands(f) + 1);
            sil_selcands(1:fixed_cands(g)) = ...
                lin_path(S - fixed_cands(f) + 2:end);
        end
    end
else
    sil_selcands = sil_circpreimage(sil_pts, ...
        sil_normals, init_sigma_sil, ...
        init_gamma, init_sigma_norm, ...
        project.cand_limits * transformed, ...
        project.cand_dists, cand_norms);
end
sil_triangles     = project.cand_ixs(sil_selcands);
sil_barycentric   = project.cand_uvs(sil_selcands, :);

for s = 1:S
    u = sil_barycentric(s, 1);
    v = sil_barycentric(s, 2);
    w = 1 - sum(sil_barycentric(s, :));
    [~, tripos] = max([u w v]);
    
    % This function assumes that 'tripos' is 2. If it isn't,
    % rotate round.
    if tripos == 1
        sil_barycentric(s, 1) = sil_barycentric(s, 2);
        sil_barycentric(s, 2) = w;
        sil_triangles(s) = ...
            mesh.edges(sil_triangles(s)).next.next.index_in_mesh;
    elseif tripos == 3
        sil_barycentric(s, 2) = sil_barycentric(s, 1);
        sil_barycentric(s, 1) = w;
        sil_triangles(s) = ...
            mesh.edges(sil_triangles(s)).next.index_in_mesh;
    end
end
    
% Pull all non-fixed points slightly away from vertices
ptsatverts = sum(abs(sil_barycentric), 2) == 0 & ...
    ~fixed_sil_pts;
sil_barycentric(ptsatverts, :) = 1e-2;
last_sil_barycen  = sil_barycentric;
last_sil_tris     = sil_triangles;
    
problem.x0(mvs:mve) = 1;

problem.objective   = @calc_energy_newImage;
problem.options     = optimset( ...
    'Jacobian',         'on'                 ...
    , 'PreCondBandwidth', Inf                  ...
    , 'Diagnostics',      'on'                 ...
    , 'Display',          'iter'               ...
    , 'OutputFcn',        @newiteration        ...
    , 'TolFun',           1e-5                 ...
    , 'MaxIter',          150                  ...
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

resultVec = result;

result = struct();
    
v          = resultVec(rvs:rve);
result.rotate   = [expm([    0 -v(3)  v(2) ;
    v(3)     0 -v(1) ;
    -v(2)  v(1)     0 ]) zeros(3, 1) ;
    zeros(1, 3) 1];

result.translate  = [ 1 0 0 resultVec(tvs) ;
    0 1 0 resultVec(tvs + 1) ;
    0 0 1 resultVec(tvs + 2) ;
    0 0 0 1              ];

result.scale = resultVec(sv);
result.shapevars = resultVec(mvs:mve);
shapevars = result.shapevars;

result.M = model.parameters.M;
result.parameters = model.parameters;
result.S = S;
result.T = T;
result.K = K;
result.conwt = conwt;
result.resnorm   = resnorm;
result.residual  = residual;
result.exitflag  = exitflag;
result.output    = output;

    function stop = newiteration(vars, ~, state)
        if strcmp(state, 'iter')
                sil_barycentric  = last_sil_barycen;
                sil_triangles    = last_sil_tris;
                update_sil_surface(vars(svs:sve));
                last_sil_barycen = sil_barycentric;
                last_sil_tris    = sil_triangles;
                last_sil_uvs     = vars(svs:sve);
        end
        stop = false;
    end

% Note: To avoid a lot of 2^-0.5 factors, this function calculates 2E,
% where E is the energy described in the paper.
    function [energy Jdata] = calc_energy_newImage(vars)
        energy = zeros(es(2) + 3 * P * (cM + 1), 1);
        vars(vs(2) + 1:end) = reshape(modes,size(vars(vs(2) + 1:end)));
        meansc = 0;
        
        % over all images -> one image of interest
        
        sil_barycentric  = last_sil_barycen;
        sil_triangles   = last_sil_tris;
        update_sil_surface(vars(svs:sve));
        
        % Calculate derivatives that give the relation between the
        % parameter space for sil_barycentric (with transitions around
        % extraordinary points), (u,v), and the uniform parameter space
        % associated with each triangle, (a,b).
        dadu = zeros(S, 1);
        dbdu = zeros(S, 1);
        dadv = zeros(S, 1);
        dbdv = zeros(S, 1);
        for j = 1:S
                if fixed_sil_pts(j)
                    continue
                end
                
                cpt     = sil_barycentric(j, 1) * u_basis + ...
                    sil_barycentric(j, 2) * v_basis;
                e       = mesh.edges(sil_triangles(j));
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
            
            v   = vars(rvs:rve);
            if nargout > 1
                silJ  = zeros(2 * S + 5 * T,            2     );
                varsJ = zeros(    S + 5 * T + 2 * K, 7 + cM);
                meshJ = zeros(    S + 5 * T + 2 * K, 3 * P );
                
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
                rs           = image.transform;
                dMdv_noscale = [ rs(1:3, 1:3) * dRdv(1:3, :) ;
                    rs(1:3, 1:3) * dRdv(4:6, :) ;
                    rs(1:3, 1:3) * dRdv(7:9, :) ];
                dMdv         = dMdv_noscale * vars(sv);
            end
            
            % Silhouette points
            rot         = [ expm([    0 -v(3)  v(2) ;
                v(3)     0 -v(1) ;
                -v(2)  v(1)     0 ]) zeros(3, 1) ;
                zeros(1, 3) 1 ];
            rotscale    = rot .* vars(sv);
            rotscale(4, 4) = 1;
            translate   = [ 1 0 0 vars(tvs    ) ;
                0 1 0 vars(tvs + 1) ;
                0 0 1 vars(tvs + 2) ;
                0 0 0 1                ];
            DT          = translate * ...
                image.transform * rotscale;
            noscale     = translate * image.transform * rot;
            pointvars   = reshape(modes * [1 ; vars(mvs:mve)], P, 3);
            transformed = [ pointvars ones(P, 1) ] * DT';
            transformed = transformed(:, 1:3);
            rotscaled   = [ pointvars ones(P, 1) ] * rotscale' * ...
                image.transform';
            rotscaled   = rotscaled(:, 1:3);
            rotated     = [ pointvars ones(P, 1) ] * rot' * ...
                image.transform';
            rotated     = rotated(:, 1:3);
            jT          = 0;
            for j = 0:S - 1
                % Continuity
                if j == S - 1, nextj = 1; else nextj = j + 2; end
                [dist ddda1 dddb1 ddda2 dddb2] = mesh.distbetween(...
                    mesh.edges(sil_triangles(j + 1)), ...
                    sil_barycentric(j + 1, 1), ...
                    sil_barycentric(j + 1, 2), ...
                    mesh.edges(sil_triangles(nextj)), ...
                    sil_barycentric(nextj, 1), ...
                    sil_barycentric(nextj, 2));
                energy(ces + j) = sqrt(2 * gamma) * dist;
                
                if nargout > 1
                    dddu1 = ddda1 * dadu(j + 1) + dddb1 * dbdu(j + 1);
                    dddv1 = ddda1 * dadv(j + 1) + dddb1 * dbdv(j + 1);
                    dddu2 = ddda2 * dadu(nextj) + dddb2 * dbdu(nextj);
                    dddv2 = ddda2 * dadv(nextj) + dddb2 * dbdv(nextj);
                    rt2g  = sqrt(2 * gamma);
                    
                    silJ(2 * T     + j + 1, 1) = rt2g * dddu1;
                    silJ(2 * T     + j + 1, 2) = rt2g * dddv1;
                    silJ(2 * T + S + j + 1, 1) = rt2g * dddu2;
                    silJ(2 * T + S + j + 1, 2) = rt2g * dddv2;
                end
                
                if fixed_sil_pts(j + 1)
                    continue;
                end
                
                [lt deriv sec] = mesh.limitevaluation(...
                    mesh.edges(sil_triangles(j + 1)), ...
                    sil_barycentric(j + 1, 1), ...
                    sil_barycentric(j + 1, 2));
                % X
                dx = lt * transformed(:, 1) - sil_pts(j + 1, 1);
                
                % Y
                dy = lt * transformed(:, 2) - sil_pts(j + 1, 2);
                
                energy(ses     + jT) = dx / model.parameters.sigma_sil;
                energy(ses + T + jT) = dy / model.parameters.sigma_sil;
                
                if nargout > 1
                    % X
                    xd = lt' * DT(1, 1:3);
                    xd = xd(:);
                    meshJ(       jT + 1, :) = xd / model.parameters.sigma_sil;
                    
                    % Y
                    yd = lt' * DT(2, 1:3);
                    yd = yd(:);
                    meshJ(T    + jT + 1, :) = yd / model.parameters.sigma_sil;
                end
                
                % U and V
                dxda    = [  0 1 0 ] * deriv;
                dxdb    = [ -1 0 0 ] * deriv;
                dxdu    = dxda * dadu(j + 1) + dxdb * dbdu(j + 1);
                dxdv    = dxda * dadv(j + 1) + dxdb * dbdv(j + 1);
                
                if nargout > 1
                    dxduval                      = dxdu * rotscaled;
                    dxdvval                      = dxdv * rotscaled;
                    silJ(jT + 1, 1)          = dxduval(1);
                    silJ(jT + 1, 2)          = dxdvval(1);
                    varsJ(jT + 1, 1:3)        = lt * pointvars * ...
                        dMdv([1 4 7], :);
                    varsJ(jT + 1, 4)          = 1;
                    varsJ(jT + 1, 7)          = lt * pointvars * ...
                        noscale(1, 1:3)';
                    
                    silJ(T + jT + 1, 1)   = dxduval(2);
                    silJ(T + jT + 1, 2)   = dxdvval(2);
                    varsJ(T + jT + 1, 1:3) = lt * pointvars * ...
                        dMdv([2 5 8], :);
                    varsJ(T + jT + 1, 5)   = 1;
                    varsJ(T + jT + 1, 7)   = lt * pointvars * ...
                        noscale(2, 1:3)';
                    
                    for m = 1:cM
                        varsJ(       jT + 1, 7 + m) = ...
                            meshJ(       jT + 1, :) * modes(:, m + 1);
                        varsJ(T + jT + 1, 7 + m) = ...
                            meshJ(T + jT + 1, :) * modes(:, m + 1);
                    end
                end
                
                % Normals
                dxduval    = dxdu * rotated;
                dxdvval    = dxdv * rotated;
                unormal    = cross(dxdvval, dxduval);
                normlen    = sqrt(sum(unormal .^ 2));
                normal     = unormal / normlen;
                s_norm_rep = 1 / model.parameters.sigma_norm;
                
                energy(nes +            jT) = s_norm_rep * (normal(1) - ...
                    sil_normals(j + 1, 1));
                energy(nes +     T + jT) = s_norm_rep * (normal(2) - ...
                    sil_normals(j + 1, 2));
                energy(nes + 2 * T + jT) = s_norm_rep * normal(3);
                
                if nargout > 1
                    cpt     = sil_barycentric(j + 1, 1) * u_basis + ...
                        sil_barycentric(j + 1, 2) * v_basis;
                    e       = mesh.edges(sil_triangles(j + 1));
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
                    
                    n1ej    = S + 2 * T + jT + 1;
                    n2ej    = S + 3 * T + jT + 1;
                    n3ej    = S + 4 * T + jT + 1;
                    
                    meshJ(n1ej, 1:P)             = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 1) + ...
                        dndy(1,:) * noscale(2, 1) + ...
                        dndz(1,:) * noscale(3, 1));
                    meshJ(n1ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 2) + ...
                        dndy(1,:) * noscale(2, 2) + ...
                        dndz(1,:) * noscale(3, 2));
                    meshJ(n1ej, 2 * P + 1:3 * P) = s_norm_rep * ...
                        (dndx(1,:) * noscale(1, 3) + ...
                        dndy(1,:) * noscale(2, 3) + ...
                        dndz(1,:) * noscale(3, 3));
                    
                    meshJ(n2ej, 1:P)             = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 1) + ...
                        dndy(2,:) * noscale(2, 1) + ...
                        dndz(2,:) * noscale(3, 1));
                    meshJ(n2ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 2) + ...
                        dndy(2,:) * noscale(2, 2) + ...
                        dndz(2,:) * noscale(3, 2));
                    meshJ(n2ej, 2 * P + 1:3 * P) = s_norm_rep * ...
                        (dndx(2,:) * noscale(1, 3) + ...
                        dndy(2,:) * noscale(2, 3) + ...
                        dndz(2,:) * noscale(3, 3));
                    
                    meshJ(n3ej, 1:P)             = s_norm_rep * ...
                        (dndx(3,:) * noscale(1, 1) + ...
                        dndy(3,:) * noscale(2, 1) + ...
                        dndz(3,:) * noscale(3, 1));
                    
                    meshJ(n3ej, P + 1:2 * P)     = s_norm_rep * ...
                        (dndx(3,:) * noscale(1, 2) + ...
                        dndy(3,:) * noscale(2, 2) + ...
                        dndz(3,:) * noscale(3, 2));
                    
                    meshJ(n3ej, 2 * P + 1:3 * P) = s_norm_rep * ...
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
                    
                    silJ(n1ej + S, 1) = s_norm_rep * dndu(1);
                    silJ(n1ej + S, 2) = s_norm_rep * dndv(1);
                    silJ(n2ej + S, 1) = s_norm_rep * dndu(2);
                    silJ(n2ej + S, 2) = s_norm_rep * dndv(2);
                    silJ(n3ej + S, 1) = s_norm_rep * dndu(3);
                    silJ(n3ej + S, 2) = s_norm_rep * dndv(3);
                    
                    % Changes to the matrix transform from rotation
                    dndvar = ...
                        (dndx * pointvars * dMdv_noscale([1 4 7], :) + ...
                        dndy * pointvars * dMdv_noscale([2 5 8], :) + ...
                        dndz * pointvars * dMdv_noscale([3 6 9], :));
                    
                    varsJ(n1ej, 1:3) = s_norm_rep * dndvar(1, :);
                    varsJ(n2ej, 1:3) = s_norm_rep * dndvar(2, :);
                    varsJ(n3ej, 1:3) = s_norm_rep * dndvar(3, :);
                    
                    for m = 1:cM
                        varsJ(n1ej, 7 + m) = ...
                            meshJ(n1ej, :) * modes(:, m + 1);
                        varsJ(n2ej, 7 + m) = ...
                            meshJ(n2ej, :) * modes(:, m + 1);
                        varsJ(n3ej, 7 + m) = ...
                            meshJ(n3ej, :) * modes(:, m + 1);
                    end
                end
                
                jT = jT + 1;
            end
            
            % Constraint vertices
            energy(cxes:cyee) = repmat(conwt, 2, 1) .* ...
                (constAeq * reshape(transformed(:, 1:2), 2 * P, 1) - ...
                constbeq);
            if nargout > 1
                C = repmat(conwt, 1, P) .* constAeq(1:K, 1:P);
                
                % Derivatives for equations measuring difference in X
                cxs = S + 5 * T + 1;
                cxe = S + 5 * T + K;
                % - Changes to mesh
                meshJ(cxs:cxe,         1:P)     = DT(1, 1) * C;
                meshJ(cxs:cxe,     P + 1:2 * P) = DT(1, 2) * C;
                meshJ(cxs:cxe, 2 * P + 1:3 * P) = DT(1, 3) * C;
                % - Changes to matrix transform
                %   - From translation
                varsJ(cxs:cxe, 4)               = conwt;
                %   - From rotation
                varsJ(cxs:cxe, 1:3)             = ...
                    C * pointvars * dMdv([1 4 7], :);
                %   - From scaling
                varsJ(cxs:cxe, 7)               = ...
                    C * pointvars * noscale(1, 1:3)';
                
                % Derivatives for equations measuring difference in Y
                cys = S + 5 * T +     K + 1;
                cye = S + 5 * T + 2 * K;
                % - Changes to mesh
                meshJ(cys:cye,         1:P)     = DT(2, 1) * C;
                meshJ(cys:cye,     P + 1:2 * P) = DT(2, 2) * C;
                meshJ(cys:cye, 2 * P + 1:3 * P) = DT(2, 3) * C;
                % - Changes to matrix transform
                %   - From translation
                varsJ(cys:cye, 5)               = conwt;
                %   - From rotation
                varsJ(cys:cye, 1:3)             = ...
                    C * pointvars * dMdv([2 5 8], :);
                %   - From scaling
                varsJ(cys:cye, 7)               = ...
                    C * pointvars * noscale(2, 1:3)';
                
                for m = 1:cM
                    varsJ(cxs:cye, 7 + m) = ...
                        meshJ(cxs:cye, :) * modes(:, m + 1);
                end
            end
            
            % Shape coefficient regularisation
            energy(mes:mee) = 2 * beta * vars(mvs:mve);
            meansc                = meansc + vars(sv);
            
            if nargout > 1
                assert(max(S) > 2 * max(K));
                Jdata = zeros(2 * max(S) + 5 * max(T) + 1, ...
                    (9 + cM + 3 * P));
                Jdata(end, 1:length(vars)) = vars';
                
                
                bs = 1;
                be = (9 + cM + 3 * P);
                Jdata(1:2 * S + 5 * T, bs:bs + 1) = silJ;
                Jdata(1:    S + 5 * T + 2 * K, ...
                    bs + 2:bs + 8 + cM) = varsJ;
                Jdata(1:    S + 5 * T + 2 * K, ...
                    bs + 9 + cM:be)     = meshJ;
                
                
                if ~jacobmult_on
                    JdataT = jacobmult(Jdata, eye(length(energy)), -1);
                    Jdata  = jacobmult(Jdata, eye(length(vars)),    1);
                    assert(~any(any(abs(JdataT' - Jdata) > 1e-10)));
                end
            end
    end

    function update_sil_surface(new_uvs)
        if all(last_sil_uvs == new_uvs)
            return
        end
        
        jT = 1;
        
        for j = 1:S
            if fixed_sil_pts(j)
                continue;
            else
                diff    =             new_uvs([jT T + jT]) - ...
                    last_sil_uvs([jT T + jT]);
                diff    = complex(diff(1), diff(2));
                jT = jT + 1;
                
                if diff == 0
                    continue;
                end
            end
            
            e           = mesh.edges(sil_triangles(j));
            last_sil_pt = sil_barycentric(j, 1) * u_basis + ...
                sil_barycentric(j, 2) * v_basis;
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
                        u = sil_barycentric(j, 1);
                        v = sil_barycentric(j, 2);
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
            sil_barycentric(j, :) = ...
                [ real(moveto) imag(moveto) ] / change_basis';
            w = 1 - sum(sil_barycentric(j, :));
            [~, tripos] = max([sil_barycentric(j, :) w]);
            if tripos == 1
                sil_barycentric(j, 1) = ...
                    sil_barycentric(j, 2);
                sil_barycentric(j, 2) = w;
                e = e.next.next;
            elseif tripos == 2
                sil_barycentric(j, 2) = ...
                    sil_barycentric(j, 1);
                sil_barycentric(j, 1) = w;
                e = e.next;
            end
            sil_triangles(j) = e.index_in_mesh;
        end
    end

    function W = jacobmult(Jdata, X, flag)
        vars = Jdata(end, 1:2 * T + (7 + cM) + 3 * P * (cM + 1));
        
        meansc = 0;
        bs       = 1;
        be       = (9 + cM + 3 * P);
        
        silJ  = Jdata(1:2 * S + 5 * T, bs:bs + 1);
        varsJ = Jdata(1:    S + 5 * T + 2 * K, ...
            bs + 2:bs + 8 + cM);
        meshJ = Jdata(1:    S + 5 * T + 2 * K, ...
            bs + 9 + cM:be);
        
        meansc   = meansc + vars(sv);
        
        if flag == 0
            flag = 1;
            both = 1;
        else
            both = 0;
        end
        
        if flag > 0
            Y = zeros(es(2) + 3 * P * (cM + 1), size(X, 2));
            
            esi = 1;
            esn = es(2) - cM;
            Y(esi:esn, :) = varsJ * X(rvs:mve, :);
            Y(esi:esn, :) = Y(esi:esn, :) + ...
                meshJ * X(vs(2) + 1:vs(2) + 3 * P, :);
            for m = 1:cM
                Y(esi:esn, :) = Y(esi:esn, :) + meshJ * ...
                    vars(mvs + m - 1) * X(vs(2) + 1 + ...
                    m * 3 * P:vs(2) + (m + 1) * 3 * P, :);
            end
            
            Y(ses:see - T, :) = Y(ses:see - T, :) ...
                + diag(silJ(1:T, 1)) * X(svs:sve - T, :);
            Y(ses:see - T, :) = Y(ses:see - T, :) ...
                + diag(silJ(1:T, 2)) * X(svs + T:sve, :);
            Y(ses + T:see, :) = Y(ses + T:see, :) ...
                + diag(silJ(T + 1:2 * T, 1)) ...
                * X(svs:sve - T, :);
            Y(ses + T:see, :) = Y(ses + T:see, :)...
                + diag(silJ(T + 1:2 * T, 2)) ...
                * X(svs + T:sve, :);
            
            D = diag(silJ(2 * T + 1:2 * T + S, 1)) + ...
                diag(silJ(2 * T + S + 1:...
                2 * T + 2 * S - 1, 1), 1);
            D(S, 1) = silJ(2 * T + 2 * S, 1);
            D(:, fixed_sil_pts) = [];
            Y(ces:cee, :) = Y(ces:cee, :) ...
                +  D * X(svs:sve - T, :);
            D = diag(silJ(2 * T + 1:2 * T + S, 2)) + ...
                diag(silJ(2 * T + S + 1:...
                2 * T + 2 * S - 1, 2), 1);
            D(S, 1) = silJ(2 * T + 2 * S, 2);
            D(:, fixed_sil_pts) = [];
            Y(ces:cee, :) = Y(ces:cee, :) ...
                + D * X(svs + T:sve, :);
            
            Y(nes:nee - 2 * T, :) = Y(nes:nee - 2 * T, :) ...
                + diag(silJ(2 * T + 2 * S + 1:3 * T + 2 * S, 1)) ...
                * X(svs:sve - T, :);
            Y(nes:nee - 2 * T, :) = Y(nes:nee - 2 * T, :) ...
                + diag(silJ(2 * T + 2 * S + 1:3 * T + 2 * S, 2)) ...
                * X(svs + T:sve, :);
            Y(nes + T:nee - T, :) = Y(nes + T:nee - T, :) ...
                + diag(silJ(3 * T + 2 * S + 1:4 * T + 2 * S, 1)) ...
                * X(svs:sve - T, :);
            Y(nes + T:nee - T, :) = Y(nes + T:nee - T, :)...
                + diag(silJ(3 * T + 2 * S + 1:4 * T + 2 * S, 2)) ...
                * X(svs + T:sve, :);
            Y(nes + 2 * T:nee, :) = Y(nes + 2 * T:nee, :) ...
                + diag(silJ(4 * T + 2 * S + 1:5 * T + 2 * S, 1)) ...
                * X(svs:sve - T, :);
            Y(nes + 2 * T:nee, :) = Y(nes + 2 * T:nee, :)...
                + diag(silJ(4 * T + 2 * S + 1:5 * T + 2 * S, 2)) ...
                * X(svs + T:sve, :);
            
            Y(es(2) + 1:es(2) + 3 * P, :) = ...
                Y(es(2) + 1:es(2) + 3 * P, :) ...
                + (model.parameters.xi_0) * tplatesqrt_3 * ...
                vars(vs(2) + 1:vs(2) + 3 * P)' * X(sv, :);
            for m = 1:cM
                ms = es(2) + 1 +       m * 3 * P;
                me = es(2) +     (m + 1) * 3 * P;
                Y(ms:me, :) = Y(ms:me, :) + tplatesqrt_3 * ...
                    (model.parameters.xi_def) * vars(vs(2) + 1 + m * 3 * P:...
                    vs(2) + (m + 1) * 3 * P)' ...
                    * X(sv, :);
            end
            
            Y(mes:mee, :) = 2 * beta * X(mvs:mve, :);
            
            Y(es(2) + 1:es(2) + 3 * P, :) = ...
                Y(es(2) + 1:es(2) + 3 * P, :) + ...
                model.parameters.xi_0 * tplatesqrt_3 * meansc * ...
                X(vs(2) + 1:vs(2) + 3 * P, :);
            for m = 1:cM
                ms = es(2) + 1 +       m * 3 * P;
                me = es(2) +     (m + 1) * 3 * P;
                
                Y(ms:me, :) = Y(ms:me, :) + model.parameters.xi_def * tplatesqrt_3 * ...
                    meansc * X(vs(2) + 1 + m * 3 * P:...
                    vs(2) + (m + 1) * 3 * P, :);
            end
            
            if both
                flag = -1;
            end
        else
            Y = X;
        end
        
        if flag < 0
            W = zeros(length(vars), size(Y, 2));
            
            esi = 1;
            esn = es(2) - cM;
            W(rvs:mve, :) = varsJ' * Y(esi:esn, :);
            
            ms  = vs(2) + 1;
            me  = vs(2) + 3 * P;
            W(ms:me, :) = W(ms:me, :) + meshJ' * Y(esi:esn, :);
            
            for m = 1:cM
                ms  = vs(2) + 1 + m * 3 * P;
                me  = vs(2) + (m + 1) * 3 * P;
                
                W(ms:me, :) = W(ms:me, :) + ...
                    vars(mvs + m - 1) * meshJ' * Y(esi:esn, :);
            end
            
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                + diag(silJ(1:T, 1)) * Y(ses:see - T, :);
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + diag(silJ(1:T, 2)) * Y(ses:see - T, :);
            
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                + diag(silJ(T + 1:2 * T, 1)) ...
                * Y(ses + T:see, :);
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + diag(silJ(T + 1:2 * T, 2)) ...
                * Y(ses + T:see, :);
            
            D = diag(silJ(2 * T + 1:2 * T + S, 1)) + ...
                diag(silJ(2 * T + S + 1:...
                2 * T + 2 * S - 1, 1), -1);
            D(1, S) = silJ(2 * T + 2 * S, 1);
            D(fixed_sil_pts, :) = [];
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                +  D * Y(ces:cee, :);
            
            D = diag(silJ(2 * T + 1:2 * T + S, 2)) + ...
                diag(silJ(2 * T + S + 1:...
                2 * T + 2 * S - 1, 2), -1);
            D(1, S) = silJ(2 * T + 2 * S, 2);
            D(fixed_sil_pts, :) = [];
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + D * Y(ces:cee, :);
            
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                + diag(silJ(2 * T + 2 * S + 1:3 * T + 2 * S, 1)) ...
                * Y(nes:nee - 2 * T, :);
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + diag(silJ(2 * T + 2 * S + 1:3 * T + 2 * S, 2)) ...
                * Y(nes:nee - 2 * T, :);
            
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                + diag(silJ(3 * T + 2 * S + 1:4 * T + 2 * S, 1)) ...
                * Y(nes + T:nee - T, :);
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + diag(silJ(3 * T + 2 * S + 1:4 * T + 2 * S, 2)) ...
                * Y(nes + T:nee - T, :);
            
            W(svs:sve - T, :) = W(svs:sve - T, :) ...
                + diag(silJ(4 * T + 2 * S + 1:5 * T + 2 * S, 1)) ...
                * Y(nes + 2 * T:nee, :);
            W(svs + T:sve, :) = W(svs + T:sve, :) ...
                + diag(silJ(4 * T + 2 * S + 1:5 * T + 2 * S, 2)) ...
                * Y(nes + 2 * T:nee, :);
            
            W(sv, :) = W(sv, :) + (model.parameters.xi_0) * ...
                vars(vs(2) + 1:vs(2) + 3 * P) * ...
                tplatesqrt_3' * Y(es(2) + 1:es(2) + 3 * P, :);
            
            for m = 1:cM
                ms = es(2) + 1 +       m * 3 * P;
                me = es(2) +     (m + 1) * 3 * P;
                
                W(sv, :) = W(sv, :) + (model.parameters.xi_def) * ...
                    vars(vs(2) + 1 + m * 3 * P:...
                    vs(2) + (m + 1) * 3 * P) * ...
                    tplatesqrt_3' * Y(ms:me, :);
            end
            
            W(mvs:mve, :) = ...
                W(mvs:mve, :) + 2 * beta * Y(mes:mee, :);
            
            W(vs(2) + 1:vs(2) + 3 * P, :) = ...
                W(vs(2) + 1:vs(2) + 3 * P, :) + ...
                model.parameters.xi_0 * meansc * tplatesqrt_3 * ...
                Y(es(2) + 1:es(2) + 3 * P, :);
            for m = 1:cM
                ms  = es(2) + 1 +      m  * 3 * P;
                me  = es(2) +     (m + 1) * 3 * P;
                sms = vs(2) + 1 +      m  * 3 * P;
                sme = vs(2) +     (m + 1) * 3 * P;
                
                W(sms:sme, :) = W(sms:sme, :) + ...
                    model.parameters.xi_def * tplatesqrt_3 * meansc * Y(ms:me, :);
            end
        else
            W = Y;
        end
        
        assert(~any(any(isnan(W))));
    end

end
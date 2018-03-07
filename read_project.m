function project = read_project(file)

% READ_PROJECT  Reads Forms project from a file

[project.ply project.images]     = read_fpj(file);
[project.vertices project.faces] = trimesh_fromply(project.ply);
project.mesh                     = meshtri(project.vertices, project.faces);
project.file                     = file;

% Required in order to init_candidates
project.mesh.allocNewVerts();
[project.cand_ixs ...
 project.cand_uvs ...
 project.cand_limits ...
 project.cand_dists ...
 project.cand_derivs]            = init_candidates(project.mesh);

end
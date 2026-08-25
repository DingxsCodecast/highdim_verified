%% Examples for both benchmark modes
generator_directory=fileparts(mfilename('fullpath'));
addpath(generator_directory);

structural=generate_highdim_dblp(20,'Q',3,'YDimension',10, ...
    'Seed',20260820,'Orientation','QR','Validate',true);
fprintf('Structural: %s\n',structural.validation.message);
fprintf('  x dimension %d, rays %.0f, degree %d, optimum %.8g\n', ...
    structural.metadata.n,structural.metadata.number_apex_neighbors, ...
    structural.metadata.cross_section_vertex_degree,structural.minimum);

lifted=generate_lifted_highdim_dblp(8,'Q',3,'Eta',0.5,'Alpha',-0.5, ...
    'Seed',20260820,'Orientation','HH','Validate',true);
fprintf('Lifted: %s\n',lifted.validation.message);
fprintf('  x dimension %d, PGM %.8g, incumbent %.8g, global %.8g\n', ...
    lifted.metadata.x_dimension,lifted.pgm_value, ...
    lifted.incumbent_value,lifted.minimum);

tri=reference_triangulation(8,3,'CoverageSamples',500,'Seed',20260820);
assert(tri.passed,tri.message);
cut=oracle_polar_cut(lifted,tri,1,'Coordinates','transformed');
assert(cut.passed);
fprintf('  reference cells %d, maximum ray-depth error %.3e\n', ...
    tri.number_cells,cut.max_relative_depth_error);

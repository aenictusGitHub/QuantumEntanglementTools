function run_symmetric_extensions_oracle()
% RUN_SYMMETRIC_EXTENSIONS_ORACLE Generate solver-free QETLAB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
% Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
required_files = { ...
    'SymmetricExtension.m', ...
    'SymmetricInnerExtension.m', ...
    'RandomPPTState.m', ...
    fullfile('helpers', 'jacobi_poly.m')};
for index = 1:length(required_files)
    if ~exist(fullfile(qetlab_path, required_files{index}), 'file')
        error('QuantumEntanglementTools:OracleConfiguration', ...
            'The selected QETLAB checkout is missing %s.', required_files{index});
    end
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

mixed = eye(4) / 4;
bell_vector = [1; 0; 0; 1] / sqrt(2);
bell = bell_vector * bell_vector';

extension_fixtures = struct( ...
    'name', {}, ...
    'matrix_name', {}, ...
    'order', {}, ...
    'dimensions', {}, ...
    'ppt', {}, ...
    'bosonic', {}, ...
    'expected', {});
extension_fixtures(end + 1) = extension_fixture( ...
    'one_copy_mixed', 'mixed', 1, [2, 2], 0, 0, ...
    SymmetricExtension(mixed, 1, [2, 2], 0, 0));
extension_fixtures(end + 1) = extension_fixture( ...
    'two_qubit_mixed_analytic', 'mixed', 2, [2, 2], 0, 0, ...
    SymmetricExtension(mixed, 2, [2, 2], 0, 0));
extension_fixtures(end + 1) = extension_fixture( ...
    'two_qubit_bell_analytic', 'bell', 2, [2, 2], 0, 0, ...
    SymmetricExtension(bell, 2, [2, 2], 0, 0));
extension_fixtures(end + 1) = extension_fixture( ...
    'low_dimension_ppt_mixed', 'mixed', 3, [2, 2], 1, 1, ...
    SymmetricExtension(mixed, 3, [2, 2], 1, 1));
extension_fixtures(end + 1) = extension_fixture( ...
    'low_dimension_ppt_bell', 'bell', 3, [2, 2], 1, 0, ...
    SymmetricExtension(bell, 3, [2, 2], 1, 0));

jacobi_fixtures = struct('name', {}, 'alpha', {}, 'beta', {}, 'degree', {}, ...
    'coefficients', {});
jacobi_fixtures(end + 1) = jacobi_fixture( ...
    'legendre_degree_two', 0, 0, 2, jacobi_poly(0, 0, 2));
jacobi_fixtures(end + 1) = jacobi_fixture( ...
    'asymmetric_degree_three', 1, 0, 3, jacobi_poly(1, 0, 3));
jacobi_fixtures(end + 1) = jacobi_fixture( ...
    'inner_ppt_order_three', 0, 1, 2, jacobi_poly(0, 1, 2));

rng(7301, 'twister');
random_state = RandomPPTState([2, 3], 6, 1e-12, 1);
random_transpose = PartialTranspose(random_state, 2, [2, 3]);
random_fixture = struct();
random_fixture.name = 'full_rank_shifted_seeded_properties';
random_fixture.seed = 7301;
random_fixture.dimensions = [2, 3];
random_fixture.requested_rank = 6;
random_fixture.trace = real(trace(random_state));
random_fixture.hermiticity_residual = norm(random_state - random_state', 'fro');
random_fixture.minimum_eigenvalue = min(real(eig(random_state)));
random_fixture.minimum_partial_transpose_eigenvalue = ...
    min(real(eig(random_transpose)));
random_fixture.numerical_rank = rank(random_state, 1e-10);
random_fixture.partial_transpose_numerical_rank = rank(random_transpose, 1e-10);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP3-symmetric-extensions-random-ppt';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.source_free_fixture = true;
metadata.symmetric_extension_source_sha256 = ...
    'ca8e1aaf9b3766a3e29bc1375eecf473cdb77bc71159a033411c77eacb6e7e59';
metadata.symmetric_inner_extension_source_sha256 = ...
    'f73e38499197fb18c1ff9cff554fd77e13ab7ccfbe48590a25437cd70d7f39da';
metadata.random_ppt_state_source_sha256 = ...
    'e5208eac2c04dac95051d8c5850bb6fdaaca6430117416d0675b2c9a26ef04fd';
metadata.jacobi_poly_source_sha256 = ...
    '3001148f6fb136306a23a48ce6176b2b2ed65b278e4f1f998cb180486124d513';
metadata.extension_fixture_count = length(extension_fixtures);
metadata.jacobi_fixture_count = length(jacobi_fixtures);
metadata.random_fixture_count = 1;
metadata.solver_scope = ...
    ['Only QETLAB solver-free analytic branches are executed. CVX-backed ', ...
     'outer and inner models are validated independently with Hypatia and SCS.'];
metadata.random_scope = ...
    ['Random entries are not compared across engines. The artifact records ', ...
     'seeded QETLAB full-rank invariants and Julia constructions are checked ', ...
     'against the same mathematical PSD, trace, PPT, and rank properties.'];
metadata.inner_warning = ...
    ['A negative SymmetricInnerExtension dual object is not automatically ', ...
     'an entanglement witness.'];
metadata.native_low_rank_deviation = ...
    ['The Julia low-rank route uses a bounded random separable mixture ', ...
     'instead of the pinned potentially unbounded pseudoinverse iteration.'];
if exist('OCTAVE_VERSION', 'builtin')
    metadata.engine = 'Octave';
    metadata.engine_version = OCTAVE_VERSION;
    metadata.matlab_compatible_oracle = false;
else
    metadata.engine = 'MATLAB';
    metadata.engine_version = version;
    metadata.matlab_compatible_oracle = true;
end

payload = struct();
payload.metadata = metadata;
payload.extension_fixtures = extension_fixtures;
payload.jacobi_fixtures = jacobi_fixtures;
payload.random_fixture = random_fixture;
encoded = jsonencode(payload);
output_directory = fileparts(output_path);
if ~isempty(output_directory) && ~exist(output_directory, 'dir')
    mkdir(output_directory);
end
file_id = fopen(output_path, 'w');
if file_id < 0
    error('QuantumEntanglementTools:OracleOutput', ...
        'Could not open the oracle output path for writing.');
end
cleanup_file = onCleanup(@() fclose(file_id));
fwrite(file_id, encoded, 'char');
fwrite(file_id, sprintf('\n'), 'char');
fprintf(['Wrote %d symmetric-extension, %d Jacobi, and 1 random-PPT ', ...
    'fixture to %s using %s %s.\n'], ...
    length(extension_fixtures), length(jacobi_fixtures), output_path, ...
    metadata.engine, metadata.engine_version);
end

function fixture = extension_fixture( ...
    name, matrix_name, order, dimensions, ppt, bosonic, expected)
fixture = struct();
fixture.name = name;
fixture.matrix_name = matrix_name;
fixture.order = order;
fixture.dimensions = dimensions;
fixture.ppt = ppt;
fixture.bosonic = bosonic;
fixture.expected = expected;
end

function fixture = jacobi_fixture(name, alpha, beta, degree, coefficients)
fixture = struct();
fixture.name = name;
fixture.alpha = alpha;
fixture.beta = beta;
fixture.degree = degree;
fixture.coefficients = reshape(coefficients, 1, []);
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

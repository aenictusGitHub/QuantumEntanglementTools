function run_induced_schatten_oracle()
% RUN_INDUCED_SCHATTEN_ORACLE Generate deterministic QETLAB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2022 Nathaniel Johnston and Mateus Araujo,
% BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
if ~exist(fullfile(qetlab_path, 'InducedSchattenNorm.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain InducedSchattenNorm.m.');
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

fixtures = struct( ...
    'name', {}, ...
    'dims', {}, ...
    'real', {}, ...
    'imaginary', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

identity_map = {eye(2)};
diagonal_map = {diag([1, 2])};
amplitude_map = {diag([1, sqrt(0.7)]); [0, sqrt(0.3); 0, 0]};
seed_oracle(1101);
[identity_value, identity_witness] = InducedSchattenNorm( ...
    identity_map, 3, 2);
seed_oracle(1102);
[diagonal_value, diagonal_witness] = InducedSchattenNorm( ...
    diagonal_map, 3, 2);
seed_oracle(1103);
[amplitude_value, amplitude_witness] = InducedSchattenNorm( ...
    amplitude_map, 4, 2);
[amplitude_exact, amplitude_exact_witness] = InducedSchattenNorm( ...
    amplitude_map, 2, 2);

fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_identity_3to2_value', identity_value);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_identity_3to2_witness', identity_witness);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_diagonal_3to2_value', diagonal_value);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_diagonal_3to2_witness', diagonal_witness);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_amplitude_4to2_value', amplitude_value);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_amplitude_4to2_witness', amplitude_witness);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_amplitude_2to2_value', amplitude_exact);
fixtures(end + 1) = oracle_fixture( ...
    'induced_schatten_amplitude_2to2_witness', amplitude_exact_witness);

tol_x0_error_identifier = '';
try
    InducedSchattenNorm( ...
        identity_map, 3, 2, [2, 2], 1e-8, eye(2));
catch tol_x0_error
    tol_x0_error_identifier = tol_x0_error.identifier;
end
if isempty(tol_x0_error_identifier)
    error('QuantumEntanglementTools:OracleUnexpectedSuccess', ...
        'The reviewed pinned TOL/X0 forwarding defect unexpectedly disappeared.');
end

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-induced-schatten';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.fixed_start = false;
metadata.seeded_random_start = true;
metadata.tol_x0_error_identifier = tol_x0_error_identifier;
metadata.evidence_scope = ...
    ['Seeded numeric branches only. The pinned optional TOL/X0 route is ', ...
     'unreachable because it forwards too many arguments to ', ...
     'superoperator_dims. Iterative agreement is lower-bound evidence and ', ...
     'does not certify a global induced norm.'];
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
payload.fixtures = fixtures;
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
fprintf('Wrote %d InducedSchattenNorm fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value)
    fixture = struct();
    fixture.name = name;
    fixture.dims = size(value);
    fixture.real = reshape(full(real(value)), 1, []);
    fixture.imaginary = reshape(full(imag(value)), 1, []);
    fixture.comparison = 'normwise';
    fixture.atol = 2e-10;
    fixture.rtol = 2e-10;
end

function seed_oracle(value)
if exist('OCTAVE_VERSION', 'builtin')
    randn('seed', value);
else
    rng(value, 'twister');
end
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

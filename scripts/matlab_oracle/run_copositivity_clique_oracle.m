function run_copositivity_clique_oracle()
% RUN_COPOSITIVITY_CLIQUE_ORACLE Generate pinned source-free fixtures.
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
required_files = {'IsCopositive.m', 'CliqueNumber.m', ...
    'CopositivePolynomial.m', 'PolynomialOptimize.m'};
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

horn = [
     1 -1  1  1 -1;
    -1  1 -1  1  1;
     1 -1  1 -1  1;
     1  1 -1  1 -1;
    -1  1  1 -1  1
];

copositivity_fixtures = struct( ...
    'name', {}, ...
    'dimension', {}, ...
    'matrix', {}, ...
    'expected', {});
copositivity_fixtures(end + 1) = copositivity_fixture( ...
    'positive_diagonal', [1 0; 0 2], 11);
copositivity_fixtures(end + 1) = copositivity_fixture( ...
    'negative_diagonal', [-1 0; 0 2], 12);
copositivity_fixtures(end + 1) = copositivity_fixture( ...
    'fixed_threshold_negative', [-5e-10 0; 0 1], 13);
copositivity_fixtures(end + 1) = copositivity_fixture( ...
    'negative_pair_ray', [1 -2; -2 1], 14);
copositivity_fixtures(end + 1) = copositivity_fixture( ...
    'horn_boundary', horn, 15);

cycle_five = [
    0 1 0 0 1;
    1 0 1 0 0;
    0 1 0 1 0;
    0 0 1 0 1;
    1 0 0 1 0
];
triangle_isolate = [
    0 1 1 0;
    1 0 1 0;
    1 1 0 0;
    0 0 0 0
];
path_four = [
    0 1 0 0;
    1 0 1 0;
    0 1 0 1;
    0 0 1 0
];

clique_fixtures = struct( ...
    'name', {}, ...
    'dimension', {}, ...
    'adjacency', {}, ...
    'expected_upper', {}, ...
    'expected_lower', {});
clique_fixtures(end + 1) = clique_fixture( ...
    'edgeless_three', zeros(3), 21);
clique_fixtures(end + 1) = clique_fixture( ...
    'complete_three', ones(3) - eye(3), 22);
clique_fixtures(end + 1) = clique_fixture( ...
    'cycle_five', cycle_five, 23);
clique_fixtures(end + 1) = clique_fixture( ...
    'triangle_plus_isolate', triangle_isolate, 24);
clique_fixtures(end + 1) = clique_fixture( ...
    'path_four', path_four, 25);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP5-copositivity-clique';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.copositivity_fixture_count = length(copositivity_fixtures);
metadata.clique_fixture_count = length(clique_fixtures);
metadata.source_free_fixture = true;
metadata.is_copositive_sha256 = ...
    '6637ea6c9d4e2f7647520087cb16b2608d5e32587bc3f87e62340ca102a93d22';
metadata.clique_number_sha256 = ...
    'dc99fa9c91a6e2be985881a3598013347acf63b42fb59303b48d8e0aac6a1eff';
metadata.copositive_polynomial_sha256 = ...
    '4814f9bfbfe4379b725035c1a6c9d60652216278b20e74bdc745c148ca0544bb';
metadata.polynomial_optimize_sha256 = ...
    'b70a02cc12c9b341a0aaf11e22a4e0382d3e5c18f9b345ff7e08791ad59a24d7';
metadata.fixed_threshold_defect = ...
    ['Pinned IsCopositive reports true for diag(-5e-10,1) because it ', ...
     'accepts hierarchy lower bounds down to -1e-9. The Julia API reports ', ...
     'this tolerance boundary as unknown.'];
metadata.evidence_scope = ...
    ['The no-SDP QETLAB path is used because CVX is unavailable. ', ...
     'Copositivity outputs are stable tri-state observations, but QETLAB ', ...
     'chooses its sample count from elapsed wall time and its global RNG. ', ...
     'The Julia comparison therefore uses exact analytic witnesses and ', ...
     'fixed explicit sample counts rather than comparing random inner values. ', ...
     'Clique fixtures compare the final integer upper/lower outputs.'];
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
payload.copositivity_fixtures = copositivity_fixtures;
payload.clique_fixtures = clique_fixtures;
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
fprintf('Wrote %d copositivity and %d clique fixtures to %s using %s %s.\n', ...
    length(copositivity_fixtures), length(clique_fixtures), output_path, ...
    metadata.engine, metadata.engine_version);
end

function fixture = copositivity_fixture(name, matrix, seed)
seed_rng(seed);
fixture = struct();
fixture.name = name;
fixture.dimension = size(matrix, 1);
fixture.matrix = reshape(full(matrix), 1, []);
fixture.expected = IsCopositive(matrix, 0, 'nosdp');
end

function fixture = clique_fixture(name, adjacency, seed)
seed_rng(seed);
[upper, lower] = CliqueNumber(adjacency, 0, 'nosdp');
fixture = struct();
fixture.name = name;
fixture.dimension = size(adjacency, 1);
fixture.adjacency = reshape(full(adjacency), 1, []);
fixture.expected_upper = upper;
fixture.expected_lower = lower;
end

function seed_rng(seed)
try
    rng(seed, 'twister');
catch
    rand('seed', seed);
    randn('seed', seed);
end
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

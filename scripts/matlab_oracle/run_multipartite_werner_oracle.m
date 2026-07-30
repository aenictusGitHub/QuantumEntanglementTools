function run_multipartite_werner_oracle()
% RUN_MULTIPARTITE_WERNER_ORACLE Generate deterministic QETLAB fixtures.
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
if ~exist(fullfile(qetlab_path, 'WernerState.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain WernerState.m.');
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

parameters = [0.05, 0.04, 0.03, 0.03, 0.02];
fixtures(end + 1) = oracle_fixture( ...
    'multipartite_werner_pinned_overwrite', WernerState(3, parameters));
fixtures(end + 1) = oracle_fixture( ...
    'werner_one_entry_vector', WernerState(3, [0.2]));

invalid_identifier = '';
try
    WernerState(2, [0.1, 0.2]);
catch invalid_error
    invalid_identifier = invalid_error.identifier;
end
if isempty(invalid_identifier)
    error('QuantumEntanglementTools:OracleUnexpectedSuccess', ...
        'The reviewed invalid multipartite parameter length unexpectedly succeeded.');
end

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-multipartite-werner';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.parameters = parameters;
metadata.pinned_defect = ...
    ['The multipartite loop reinitializes rho at each permutation and ', ...
     'retains only the final alpha/permutation term.'];
metadata.invalid_identifier = invalid_identifier;
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
fprintf('Wrote %d WernerState fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value)
    fixture = struct();
    fixture.name = name;
    fixture.dims = size(value);
    fixture.real = reshape(full(real(value)), 1, []);
    fixture.imaginary = reshape(full(imag(value)), 1, []);
    fixture.comparison = 'normwise';
    fixture.atol = 3e-14;
    fixture.rtol = 3e-14;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

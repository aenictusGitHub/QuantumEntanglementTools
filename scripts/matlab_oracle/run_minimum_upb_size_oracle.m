function run_minimum_upb_size_oracle()
% RUN_MINIMUM_UPB_SIZE_ORACLE Generate deterministic QETLAB fixtures.
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
if ~exist(fullfile(qetlab_path, 'MinUPBSize.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain MinUPBSize.m.');
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

fixtures(end + 1) = oracle_fixture('minimum_upb_2x3', MinUPBSize([2, 3], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_3x3x3', MinUPBSize([3, 3, 3], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_2x3x3', MinUPBSize([2, 3, 3], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_4x6', MinUPBSize([4, 6], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_four_qubits', MinUPBSize(2 * ones(1, 4), 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_six_qubits', MinUPBSize(2 * ones(1, 6), 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_eight_qubits', MinUPBSize(2 * ones(1, 8), 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_twelve_qubits', MinUPBSize(2 * ones(1, 12), 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_2x2x3', MinUPBSize([2, 2, 3], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_2x2x5', MinUPBSize([2, 2, 5], 0));
fixtures(end + 1) = oracle_fixture('minimum_upb_2x2x9', MinUPBSize([2, 2, 9], 0));

unknown_identifier = '';
try
    MinUPBSize([2, 3, 4], 0);
catch unknown_error
    unknown_identifier = unknown_error.identifier;
end
if isempty(unknown_identifier)
    error('QuantumEntanglementTools:OracleUnexpectedSuccess', ...
        'The reviewed unknown MinUPBSize case unexpectedly returned a value.');
end

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-minimum-upb-size';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.unknown_dimensions = [2, 3, 4];
metadata.unknown_identifier = unknown_identifier;
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
fprintf('Wrote %d MinUPBSize fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value)
    fixture = struct();
    fixture.name = name;
    fixture.dims = size(value);
    fixture.real = reshape(full(real(value)), 1, []);
    fixture.imaginary = reshape(full(imag(value)), 1, []);
    fixture.comparison = 'exact';
    fixture.atol = 0;
    fixture.rtol = 0;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

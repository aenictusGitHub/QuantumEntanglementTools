function run_twirl_oracle()
% RUN_TWIRL_ORACLE Generate deterministic QETLAB Twirl fixtures.
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
if ~exist(fullfile(qetlab_path, 'Twirl.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain Twirl.m.');
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

X2 = reshape(1:16, 4, 4) + 1i * reshape([0:7, -8:-1], 4, 4);
X3 = reshape(1:64, 8, 8) + 1i * reshape([0:31, -32:-1], 8, 8);
fixtures(end + 1) = oracle_fixture('werner_p2', Twirl(X2, 'werner', 2));
fixtures(end + 1) = oracle_fixture('isotropic_p2', Twirl(X2, 'isotropic', 2));
fixtures(end + 1) = oracle_fixture('real_p2', Twirl(X2, 'real', 2));
fixtures(end + 1) = oracle_fixture('pauli_p2', Twirl(X2, 'pauli', 2));
fixtures(end + 1) = oracle_fixture('werner_p3', Twirl(X3, 'werner', 3));
fixtures(end + 1) = oracle_fixture('real_p3', Twirl(X3, 'real', 3));

invalid_type_identifier = caught_identifier(@() Twirl(X2, 'unknown', 2));
invalid_isotropic_copies_identifier = caught_identifier( ...
    @() Twirl(X3, 'isotropic', 3));
invalid_pauli_dimension_identifier = caught_identifier( ...
    @() Twirl(eye(9), 'pauli', 2));

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-twirl';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.source_sha256 = ...
    'a2ed4de3c937cb3e451d0b8e84508cdbe4b3bc23013558015d0e2b972853cd33';
metadata.invalid_type_identifier = invalid_type_identifier;
metadata.invalid_isotropic_copies_identifier = invalid_isotropic_copies_identifier;
metadata.invalid_pauli_dimension_identifier = invalid_pauli_dimension_identifier;
metadata.corrected_copy_validation = ...
    ['The pinned routine rejects isotropic/Pauli P only when P > 2; ', ...
     'the Julia API requires P = 2 before applying either bipartite formula.'];
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
fprintf('Wrote %d Twirl fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function identifier = caught_identifier(callback)
identifier = '';
try
    callback();
catch caught_error
    identifier = caught_error.identifier;
end
if isempty(identifier)
    error('QuantumEntanglementTools:OracleUnexpectedSuccess', ...
        'A reviewed invalid Twirl call did not report an identifier.');
end
end

function fixture = oracle_fixture(name, value)
fixture = struct();
fixture.name = name;
fixture.dims = size(value);
fixture.real = reshape(full(real(value)), 1, []);
fixture.imaginary = reshape(full(imag(value)), 1, []);
fixture.comparison = 'normwise';
fixture.atol = 8e-13;
fixture.rtol = 8e-13;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

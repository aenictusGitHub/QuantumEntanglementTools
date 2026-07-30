function run_entangling_gate_oracle()
% RUN_ENTANGLING_GATE_ORACLE Generate deterministic QETLAB fixtures.
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
if ~exist(fullfile(qetlab_path, 'IsEntanglingGate.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain IsEntanglingGate.m.');
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

identity_gate = eye(4);
swap_gate = SwapOperator(2, 1);
cnot_gate = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0];
controlled_z = diag([1, 1, 1, -1]);

identity_flag = IsEntanglingGate(identity_gate);
swap_flag = IsEntanglingGate(swap_gate);
[cnot_flag, cnot_witness] = IsEntanglingGate(cnot_gate);
[controlled_z_flag, controlled_z_candidate] = ...
    IsEntanglingGate(controlled_z);

fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_identity_flag', identity_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_swap_flag', swap_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_cnot_flag', cnot_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_cnot_witness', cnot_witness);
fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_controlled_z_flag', controlled_z_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_entangling_gate_controlled_z_last_candidate', ...
    controlled_z_candidate);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-entangling-gate';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.witness_review = ...
    ['The pinned controlled-Z flag is true, but its returned final ', ...
     'two-support candidate is unnormalized, has a product output, and is ', ...
     'not a valid witness.'];
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
fprintf('Wrote %d IsEntanglingGate fixtures to %s using %s %s.\n', ...
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

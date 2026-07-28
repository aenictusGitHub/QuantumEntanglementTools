function run_tier_a_oracle()
% RUN_TIER_A_ORACLE Generate deterministic, legally safe QETLAB fixtures.
%
% This development-only script is invoked by run_tier_a_oracle.sh. It reads
% QET_ORACLE_QETLAB_PATH, QET_ORACLE_OUTPUT, and QET_ORACLE_QETLAB_COMMIT
% from the environment, adds the pinned checkout without modifying it, and
% writes a JSON document. Complex arrays are serialized as separate real and
% imaginary column-major vectors together with their shape.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
% Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');

if isempty(qetlab_path)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QET_ORACLE_QETLAB_PATH is required.');
end
if isempty(output_path)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QET_ORACLE_OUTPUT is required.');
end
if isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QET_ORACLE_QETLAB_COMMIT is required.');
end
if ~exist(fullfile(qetlab_path, 'PartialTrace.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain PartialTrace.m.');
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

fixtures = struct( ...
    'name', {}, ...
    'rows', {}, ...
    'columns', {}, ...
    'real', {}, ...
    'imaginary', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

a2 = [1, 2; 3, 4];
b2 = [0, 1; 1, 0];
a3 = [0, 1, 0; 1, 0, 1; 0, 1, 0];

fixtures(end + 1) = oracle_fixture( ...
    'tensor', Tensor(a2, b2), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'tensor_sum', TensorSum([2; -1], eye(2), a2), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'kronecker_sum', KroneckerSum(a2, a3), 'exact', 0, 0);

dims = [2, 3, 2];
permutation = [3, 1, 2];
vector12 = (1:12)';
fixtures(end + 1) = oracle_fixture( ...
    'permute_systems', ...
    PermuteSystems(vector12, permutation, dims), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'permutation_operator', ...
    PermutationOperator(dims, permutation), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'swap', Swap(vector12, [1, 3], dims), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'swap_operator', SwapOperator([2, 3]), 'exact', 0, 0);

matrix12 = reshape(1:144, 12, 12);
fixtures(end + 1) = oracle_fixture( ...
    'partial_trace', ...
    PartialTrace(matrix12, [1, 3], dims), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'partial_transpose', ...
    PartialTranspose(matrix12, [1, 3], dims), ...
    'exact', 0, 0);

matrix6 = reshape(1:36, 6, 6);
fixtures(end + 1) = oracle_fixture( ...
    'realignment', Realignment(matrix6, [2, 3]), 'exact', 0, 0);

complex_matrix = reshape((1:16) + 1i * (16:-1:1), 4, 4);
fixtures(end + 1) = oracle_fixture( ...
    'complex_partial_transpose', ...
    PartialTranspose(complex_matrix, 1, [2, 2]), ...
    'exact', 0, 0);

fixtures(end + 1) = oracle_fixture( ...
    'symmetric_projection', ...
    SymmetricProjection(3, 2, 0), ...
    'normwise', 1e-14, 1e-14);
fixtures(end + 1) = oracle_fixture( ...
    'antisymmetric_projection', ...
    AntisymmetricProjection(3, 2, 0), ...
    'normwise', 1e-14, 1e-14);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);

if exist('OCTAVE_VERSION', 'builtin')
    metadata.engine = 'Octave';
    metadata.engine_version = OCTAVE_VERSION;
    metadata.matlab_compatible_oracle = false;
else
    metadata.engine = 'MATLAB';
    metadata.engine_version = version;
    metadata.matlab_compatible_oracle = true;
end

metadata.cvx_version = 'not detected';
metadata.cvx_solver = 'not detected';
if exist('cvx_version', 'file')
    try
        metadata.cvx_version = cvx_version;
    catch
        metadata.cvx_version = 'detected but version query failed';
    end
end
if exist('cvx_solver', 'file')
    try
        metadata.cvx_solver = cvx_solver;
    catch
        metadata.cvx_solver = 'detected but solver query failed';
    end
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

fprintf('Wrote %d Tier A fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value, comparison, atol, rtol)
    value = full(value);
    fixture = struct();
    fixture.name = name;
    fixture.rows = size(value, 1);
    fixture.columns = size(value, 2);
    fixture.real = reshape(real(value), 1, []);
    fixture.imaginary = reshape(imag(value), 1, []);
    fixture.comparison = comparison;
    fixture.atol = atol;
    fixture.rtol = rtol;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

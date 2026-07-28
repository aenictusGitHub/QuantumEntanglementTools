function run_tier_c_oracle()
% RUN_TIER_C_ORACLE Generate deterministic Tier C QETLAB fixtures.
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
if ~exist(fullfile(qetlab_path, 'ApplyMap.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain ApplyMap.m.');
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

depolarizing = DepolarizingChannel(3, 0.2);
dephasing = DephasingChannel(3, 0.25);
pauli = PauliChannel([0.1, 0.2, 0.3, 0.4]);
choi_positive = ChoiMap(1, 1, 0);
reduction = ReductionMap(3, 2);

fixtures(end + 1) = oracle_fixture( ...
    'depolarizing_channel', depolarizing, 'normwise', 4e-14, 4e-14);
fixtures(end + 1) = oracle_fixture( ...
    'dephasing_channel', dephasing, 'normwise', 4e-14, 4e-14);
fixtures(end + 1) = oracle_fixture( ...
    'pauli_channel', pauli, 'normwise', 4e-14, 4e-14);
fixtures(end + 1) = oracle_fixture( ...
    'choi_map', choi_positive, 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'reduction_map', reduction, 'exact', 0, 0);

input = [0.6, 0.1 + 0.2i, 0.0; ...
         0.1 - 0.2i, 0.3, 0.05i; ...
         0.0, -0.05i, 0.1];
fixtures(end + 1) = oracle_fixture( ...
    'apply_map', ApplyMap(input, depolarizing), ...
    'normwise', 4e-14, 4e-14);

product_input = kron([0.7, 0.1i; -0.1i, 0.3], ...
                     [0.4, 0.05; 0.05, 0.6]);
fixtures(end + 1) = oracle_fixture( ...
    'partial_map', PartialMap(product_input, DephasingChannel(2), 2, [2, 2]), ...
    'normwise', 4e-14, 4e-14);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'C';
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
fprintf('Wrote %d Tier C fixtures to %s using %s %s.\n', ...
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

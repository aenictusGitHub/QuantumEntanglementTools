function run_tier_b_oracle()
% RUN_TIER_B_ORACLE Generate deterministic Tier B QETLAB fixtures.
%
% Complex arrays are serialized as separate real/imaginary column-major
% vectors. The pinned checkout is added read-only and removed on exit.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014 Nathaniel Johnston and named coauthors,
% BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');

if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
if ~exist(fullfile(qetlab_path, 'Pauli.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain Pauli.m.');
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

fixtures(end + 1) = oracle_fixture( ...
    'pauli', Pauli([1, 2, 3], 0), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'generalized_pauli', GenPauli(1, 2, 3), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'gell_mann', GellMann(8), 'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'generalized_gell_mann', GenGellMann(2, 1, 4), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'fourier_matrix', FourierMatrix(4), 'normwise', 3e-14, 3e-14);

fixtures(end + 1) = oracle_fixture( ...
    'maximally_entangled', MaxEntangled(3), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'bell_state', Bell(5), 'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'ghz_state', GHZState(3, 2, [1, 2i, -3]), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'w_state', WState(4, [1, 2, 3, 4]), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'dicke_state', DickeState(5, 2), 'normwise', 3e-14, 3e-14);

fixtures(end + 1) = oracle_fixture( ...
    'isotropic_state', IsotropicState(3, 0.25), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'werner_state', WernerState(3, 0.2), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'horodecki_3x3', HorodeckiState(0.3, [3, 3]), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'horodecki_2x4', HorodeckiState(0.3, [2, 4]), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'gisin_state', GisinState(0.4, 0.7), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'breuer_state', BreuerState(4, 0.35), ...
    'normwise', 3e-14, 3e-14);
fixtures(end + 1) = oracle_fixture( ...
    'brauer_states', BrauerStates(2, 2), 'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'chessboard_state', ChessboardState(1, 2, 3, 4, 5, 6), ...
    'normwise', 3e-14, 3e-14);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'B';
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

fprintf('Wrote %d Tier B fixtures to %s using %s %s.\n', ...
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

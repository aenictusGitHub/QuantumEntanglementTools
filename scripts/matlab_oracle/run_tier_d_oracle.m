function run_tier_d_oracle()
% RUN_TIER_D_ORACLE Generate deterministic Tier D QETLAB fixtures.
%
% Only scalar measures, singular values, and boolean necessary-condition
% evidence with unambiguous native mappings are serialized. Schmidt vectors
% are deliberately excluded because singular-vector phases and degenerate
% bases are not canonical.
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
if ~exist(fullfile(qetlab_path, 'TraceNorm.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain TraceNorm.m.');
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

operator = [1 + 2i, -2, 0.5i; ...
            3, -1i, 4 - 2i; ...
            0.25, 2 + 1i, -3];
fixtures(end + 1) = oracle_fixture( ...
    'trace_norm', TraceNorm(operator), 'normwise', 5e-12, 5e-12);
fixtures(end + 1) = oracle_fixture( ...
    'schatten_norm_p3', SchattenNorm(operator, 3), ...
    'normwise', 5e-12, 5e-12);
fixtures(end + 1) = oracle_fixture( ...
    'ky_fan_norm_k2', KyFanNorm(operator, 2), ...
    'normwise', 5e-12, 5e-12);

diagonal_state = diag([0.5, 0.3, 0.2]);
fixtures(end + 1) = oracle_fixture( ...
    'purity', Purity(diagonal_state), 'normwise', 5e-14, 5e-14);
fixtures(end + 1) = oracle_fixture( ...
    'entropy_base3_alpha1', Entropy(diagonal_state, 3, 1), ...
    'normwise', 5e-13, 5e-13);

rho = [0.7, 0.1i; -0.1i, 0.3];
sigma = [0.4, 0.05; 0.05, 0.6];
fixtures(end + 1) = oracle_fixture( ...
    'fidelity_root', Fidelity(rho, sigma), ...
    'normwise', 5e-12, 5e-12);

bell = [1; 0; 0; 1] / sqrt(2);
bell_density = bell * bell';
fixtures(end + 1) = oracle_fixture( ...
    'negativity_bell', Negativity(bell_density, [2, 2]), ...
    'normwise', 5e-13, 5e-13);

schmidt_vector = [1; 2i; 3; 4i; 5; 6i];
fixtures(end + 1) = oracle_fixture( ...
    'schmidt_coefficients_2x3', ...
    SchmidtDecomposition(schmidt_vector, [2, 3]), ...
    'normwise', 5e-12, 5e-12);
fixtures(end + 1) = oracle_fixture( ...
    'schmidt_rank_2x3', ...
    SchmidtRank(schmidt_vector, [2, 3], 1e-10), ...
    'exact', 0, 0);

mixed_bell = 0.7 * bell_density + 0.3 * eye(4) / 4;
fixtures(end + 1) = oracle_fixture( ...
    'concurrence_mixed_bell', Concurrence(mixed_bell), ...
    'normwise', 5e-12, 5e-12);

fixtures(end + 1) = oracle_fixture( ...
    'is_ppt_bell_detected', IsPPT(bell_density, 2, [2, 2], 1e-10), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'is_ppt_maximally_mixed', IsPPT(eye(4) / 4, 2, [2, 2], 1e-10), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'realignment_trace_norm_bell', ...
    TraceNorm(Realignment(bell_density, [2, 2])), ...
    'normwise', 5e-12, 5e-12);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'D';
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

fprintf('Wrote %d Tier D fixtures to %s using %s %s.\n', ...
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

function run_sk_norms_oracle()
% RUN_SK_NORMS_ORACLE Generate solver-free QETLAB norm fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014-2022 Nathaniel Johnston, BSD-2-Clause.
% Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
required_files = {'kpNormDual.m', 'SkOperatorNorm.m', 'IsBlockPositive.m'};
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

dual_fixtures = struct( ...
    'name', {}, 'kind', {}, 'k', {}, 'p', {}, 'dimensions', {}, ...
    'real', {}, 'imaginary', {}, 'expected', {});
dual_fixtures(end + 1) = dual_fixture( ...
    'vector_k2_p2', 'vector', 2, 2, [3; 2; 1]);
dual_fixtures(end + 1) = dual_fixture( ...
    'matrix_k2_p1', 'matrix', 2, 1, diag([3, 2, 1]));
dual_fixtures(end + 1) = dual_fixture( ...
    'matrix_k2_pinf', 'matrix', 2, Inf, diag([3, 2, 1]));
dual_fixtures(end + 1) = dual_fixture( ...
    'vector_k1_p3', 'vector', 1, 3, [3; 2; 1]);

sk_fixtures = struct( ...
    'name', {}, 'matrix_name', {}, 'k', {}, 'dimensions', {}, ...
    'strength', {}, 'lower', {}, 'upper', {});
diagonal = diag([4, 3, 2, 1]);
[diagonal_lower, ~, diagonal_upper, ~] = ...
    SkOperatorNorm(diagonal, 2, [2, 2], 0);
sk_fixtures(end + 1) = sk_fixture( ...
    'full_k_diagonal', 'diagonal', 2, [2, 2], 0, ...
    diagonal_lower, diagonal_upper);

bell = [1; 0; 0; 1] / sqrt(2);
product = [1; 0; 0; 0];
rank_one = 2 * bell * product';
[rank_one_lower, ~, rank_one_upper, ~] = ...
    SkOperatorNorm(rank_one, 1, [2, 2], 0);
sk_fixtures(end + 1) = sk_fixture( ...
    'rank_one_bell_product', 'rank_one', 1, [2, 2], 0, ...
    rank_one_lower, rank_one_upper);

block_fixtures = struct( ...
    'name', {}, 'matrix_name', {}, 'k', {}, 'dimensions', {}, ...
    'strength', {}, 'expected', {});
block_fixtures(end + 1) = block_fixture( ...
    'positive_identity_full_k', 'identity', 2, [2, 2], 0, ...
    IsBlockPositive(eye(4), 2, [2, 2], 0));
block_fixtures(end + 1) = block_fixture( ...
    'negative_product_diagonal', 'negative_diagonal', 1, [2, 2], 0, ...
    IsBlockPositive(diag([-1, 2, 2, 2]), 1, [2, 2], 0));
boundary = 0.5 * eye(4) - bell * bell';
block_fixtures(end + 1) = block_fixture( ...
    'bell_witness_boundary_unknown', 'bell_boundary', 1, [2, 2], 0, ...
    IsBlockPositive(boundary, 1, [2, 2], 0));
block_fixtures(end + 1) = block_fixture( ...
    'negative_entangled_full_k', 'negative_bell', 2, [2, 2], 0, ...
    IsBlockPositive(-bell * bell', 2, [2, 2], 0));

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP7-sk-norms-block-positivity';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.source_free_fixture = true;
metadata.kp_norm_dual_source_sha256 = ...
    'ecedafbbb1a58977b9447006942c23f729e79b8445139733aa6fb3ee2da2730b';
metadata.sk_operator_norm_source_sha256 = ...
    '67fc32b39abb01fb1c098750b19de1ef0cadcf622d610ec9446d9d2a0d87eaac';
metadata.is_block_positive_source_sha256 = ...
    'cfb38c34eacb699e8eb5760f1d1e6257c4f0b4a574f3713f0c3d6a2472d575a1';
metadata.dual_fixture_count = length(dual_fixtures);
metadata.sk_fixture_count = length(sk_fixtures);
metadata.block_fixture_count = length(block_fixtures);
metadata.solver_scope = ...
    ['Only numeric kpNormDual and analytic solver-free SkOperatorNorm/', ...
     'IsBlockPositive branches are executed. CVX is not invoked.'];
metadata.random_scope = ...
    ['All SkOperatorNorm fixtures terminate before sk_iterate, so no random ', ...
     'stream is consumed. Julia randomized searches have separate explicit-', ...
     'RNG property tests.'];
metadata.boundary_scope = ...
    ['The Bell witness at the exact product-state boundary is expected to ', ...
     'remain inconclusive under the pinned strict tolerance comparison.'];
metadata.model_expression_scope = ...
    ['The affine kpNormDual epigraph is validated independently with ', ...
     'Hypatia and SCS because the pinned CVX expression branch is unavailable.'];
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
payload.dual_fixtures = dual_fixtures;
payload.sk_fixtures = sk_fixtures;
payload.block_fixtures = block_fixtures;
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
fprintf(['Wrote %d dual-norm, %d S(k)-norm, and %d block-positivity ', ...
    'fixtures to %s using %s %s.\n'], ...
    length(dual_fixtures), length(sk_fixtures), length(block_fixtures), ...
    output_path, metadata.engine, metadata.engine_version);
end

function fixture = dual_fixture(name, kind, k, p, value)
fixture = struct();
fixture.name = name;
fixture.kind = kind;
fixture.k = k;
if isinf(p)
    fixture.p = 'Inf';
else
    fixture.p = p;
end
fixture.dimensions = size(value);
fixture.real = reshape(full(real(value)), 1, []);
fixture.imaginary = reshape(full(imag(value)), 1, []);
fixture.expected = kpNormDual(value, k, p);
end

function fixture = sk_fixture( ...
    name, matrix_name, k, dimensions, strength, lower, upper)
fixture = struct();
fixture.name = name;
fixture.matrix_name = matrix_name;
fixture.k = k;
fixture.dimensions = dimensions;
fixture.strength = strength;
fixture.lower = lower;
fixture.upper = upper;
end

function fixture = block_fixture( ...
    name, matrix_name, k, dimensions, strength, expected)
fixture = struct();
fixture.name = name;
fixture.matrix_name = matrix_name;
fixture.k = k;
fixture.dimensions = dimensions;
fixture.strength = strength;
fixture.expected = expected;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

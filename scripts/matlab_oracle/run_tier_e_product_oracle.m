function run_tier_e_product_oracle()
% RUN_TIER_E_PRODUCT_ORACLE Generate deterministic product-analysis fixtures.
%
% Only operator Schmidt coefficients/ranks and scalar classifications are
% serialized. Singular vectors and product factors are deliberately excluded
% because their phases, scalings, and degenerate bases are not canonical.
%
% Error flags preserve two reviewed failures in the pinned
% OperatorSchmidtDecomposition implementation. A separate flag preserves its
% NaN result for zero-concurrence mixed-state entanglement of formation.
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
if ~exist(fullfile(qetlab_path, 'OperatorSchmidtDecomposition.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain OperatorSchmidtDecomposition.m.');
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

e11 = [1, 0; 0, 0];
e22 = [0, 0; 0, 1];
basis_12 = zeros(3);
basis_12(1, 2) = 1;
basis_23 = zeros(3);
basis_23(2, 3) = 1;
operator = 3 * kron(e11, basis_12) + 2 * kron(e22, basis_23);

% K=4 requests the complete coefficient vector while avoiding phase-dependent
% factors. The non-Hermitian square input also avoids the reviewed Hermitian
% postprocessing failure recorded below.
coefficients = OperatorSchmidtDecomposition(operator, [2, 3], 4);
fixtures(end + 1) = oracle_fixture( ...
    'operator_schmidt_coefficients_2x3', coefficients, ...
    'normwise', 5e-12, 5e-12);
fixtures(end + 1) = oracle_fixture( ...
    'operator_schmidt_rank_2x3', OperatorSchmidtRank(operator, [2, 3]), ...
    'exact', 0, 0);

% Although rectangular dimensions are documented, the pinned implementation's
% unconditional Hermiticity predicate subtracts X' from rectangular X.
rectangular_left = [1, 0, 0; 0, 0, 0];
rectangular_operator = kron(rectangular_left, e11);
rectangular_error = 0;
try
    OperatorSchmidtDecomposition( ...
        rectangular_operator, [2, 2; 3, 2], 4);
catch
    rectangular_error = 1;
end
fixtures(end + 1) = oracle_fixture( ...
    'operator_schmidt_rectangular_qetlab_error', rectangular_error, ...
    'exact', 0, 0);

% For a Hermitian 2-by-3 operator, the pinned repair branch linearly indexes
% the 2-by-2 dimension table and constructs two 2-dimensional local bases.
% The resulting basis has the wrong size for the 6-by-6 input.
hermitian_operator = diag(1:6);
hermitian_error = 0;
try
    OperatorSchmidtDecomposition(hermitian_operator, [2, 3], 4);
catch
    hermitian_error = 1;
end
fixtures(end + 1) = oracle_fixture( ...
    'operator_schmidt_hermitian_2x3_qetlab_error', hermitian_error, ...
    'exact', 0, 0);

product_vector = kron(kron([1; 2i], [3; -1]), [2; 1i]);
bell = [1; 0; 0; 1] / sqrt(2);
fixtures(end + 1) = oracle_fixture( ...
    'is_product_vector_three_party_product', ...
    IsProductVector(product_vector, [2, 2, 2]), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'is_product_vector_bell_nonproduct', ...
    IsProductVector(bell, [2, 2]), ...
    'exact', 0, 0);

first_operator = [1, 2; 3, 4];
second_operator = [0, 1, 2; 3, 4, 5];
third_operator = [1; 2];
product_operator = kron(kron(first_operator, second_operator), third_operator);
rank_two_operator = kron(e11, e11) + kron(e22, e22);
fixtures(end + 1) = oracle_fixture( ...
    'is_product_operator_rectangular_three_party_product', ...
    IsProductOperator(product_operator, [2, 2, 2; 2, 3, 1]), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'is_product_operator_rank_two_nonproduct', ...
    IsProductOperator(rank_two_operator, [2, 2]), ...
    'exact', 0, 0);

angle = 0.31;
rectangular_pure = [cos(angle); 0; 0; 0; sin(angle); 0];
fixtures(end + 1) = oracle_fixture( ...
    'entformation_pure_2x3', ...
    EntFormation(rectangular_pure, [2, 3]), ...
    'normwise', 5e-12, 5e-12);

bell_density = bell * bell';
mixed_bell = 0.7 * bell_density + 0.3 * eye(4) / 4;
fixtures(end + 1) = oracle_fixture( ...
    'entformation_mixed_bell', ...
    EntFormation(mixed_bell, [2, 2]), ...
    'normwise', 5e-12, 5e-12);

% C=0 makes one binary-entropy probability zero. The pinned formula evaluates
% 0*log2(0) directly and returns NaN rather than the limiting value zero.
zero_concurrence_result = EntFormation(eye(4) / 4, [2, 2]);
fixtures(end + 1) = oracle_fixture( ...
    'entformation_zero_concurrence_qetlab_nan', ...
    double(isnan(zero_concurrence_result)), ...
    'exact', 0, 0);

fixtures(end + 1) = oracle_fixture( ...
    'in_separable_ball_maximally_mixed', ...
    InSeparableBall(eye(4) / 4), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'in_separable_ball_pure_spectrum_outside', ...
    InSeparableBall([1; 0; 0; 0]), ...
    'exact', 0, 0);

% The pinned entry point divides positive inputs by their trace. The native
% API deliberately validates trace-one input instead of normalizing it.
fixtures(end + 1) = oracle_fixture( ...
    'in_separable_ball_unnormalized_qetlab_normalizes', ...
    InSeparableBall(2 * eye(4)), ...
    'exact', 0, 0);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'E-product-analysis';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.reviewed_operator_schmidt_discrepancy = ...
    ['The pinned decomposition errors for rectangular operators and for ', ...
     'Hermitian operators with unequal local dimensions'];
metadata.reviewed_entformation_discrepancy = ...
    'The pinned mixed-state formula returns NaN at zero concurrence';
metadata.product_classification_semantics = ...
    ['Only classifications with a clear numerical margin are compared; ', ...
     'native structured boundary status is not coerced to a QETLAB boolean'];
metadata.separable_ball_semantics = ...
    ['QETLAB false means outside a sufficient ball, not entangled; the ', ...
     'native API also refuses QETLAB-style implicit trace normalization'];

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

fprintf('Wrote %d Tier E product-analysis fixtures to %s using %s %s.\n', ...
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

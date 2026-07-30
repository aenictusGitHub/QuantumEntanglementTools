function run_state_discrimination_oracle()
% RUN_STATE_DISCRIMINATION_ORACLE Generate deterministic QETLAB fixtures.
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
if ~exist(fullfile(qetlab_path, 'Distinguishability.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain Distinguishability.m.');
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

ket0 = [1; 0];
ket1 = [0; 1];
ket_plus = [1; 1] / sqrt(2);
rho0 = ket0 * ket0';
rho1 = ket1 * ket1';
rho_plus = ket_plus * ket_plus';

fixtures = struct( ...
    'name', {}, ...
    'input_kind', {}, ...
    'dimension', {}, ...
    'state_count', {}, ...
    'state_dims', {}, ...
    'state_real', {}, ...
    'state_imaginary', {}, ...
    'priors', {}, ...
    'expected_value', {}, ...
    'measurement_real', {}, ...
    'measurement_imaginary', {}, ...
    'measurement_completeness_residual', {}, ...
    'measurement_positivity_violation', {}, ...
    'measurement_objective', {}, ...
    'measurement_objective_residual', {}, ...
    'atol', {}, ...
    'rtol', {});

fixtures(end + 1) = oracle_fixture( ...
    'pure_orthogonal_equal', 'pure_columns', [ket0, ket1], [0.5, 0.5]);
fixtures(end + 1) = oracle_fixture( ...
    'pure_nonorthogonal_equal', 'pure_columns', ...
    [ket0, ket_plus], [0.5, 0.5]);
fixtures(end + 1) = oracle_fixture( ...
    'pure_nonorthogonal_unequal', 'pure_columns', ...
    [ket0, ket_plus], [0.8, 0.2]);
fixtures(end + 1) = oracle_fixture( ...
    'density_mixed_equal', 'density_stack', ...
    {diag([0.75, 0.25]), diag([0.25, 0.75])}, [0.5, 0.5]);
fixtures(end + 1) = oracle_fixture( ...
    'density_identical_unequal', 'density_stack', ...
    {rho_plus, rho_plus}, [0.8, 0.2]);
fixtures(end + 1) = oracle_fixture( ...
    'pure_orthogonal_three', 'pure_columns', eye(3), ones(1, 3) / 3);
fixtures(end + 1) = oracle_fixture( ...
    'density_orthogonal_three', 'density_stack', ...
    {diag([1, 0, 0]), diag([0, 1, 0]), diag([0, 0, 1])}, ...
    ones(1, 3) / 3);
fixtures(end + 1) = oracle_fixture( ...
    'single_pure_state', 'pure_columns', ket0, 1);

% Record two reviewed unsafe upstream behaviors as discrepancy evidence.
[pure_normalized_value, ~] = Distinguishability([2 * ket0, 3 * ket1]);
[density_normalized_value, ~] = Distinguishability({2 * rho0, 3 * rho1});
[negative_prior_value, ~] = Distinguishability({rho0, rho1}, [1.1, -0.1]);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP4-state-discrimination';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.distinguishability_sha256 = ...
    '714b3fe26124fa08520f19ed363f7cf4c856d752f6ce8b746dfc282715f111f6';
metadata.normalize_cols_sha256 = ...
    '9c1f28f57a263233e1f947a27e7ffde896a50e4c2b5294f7aff6795b68c2310d';
metadata.pure_silent_normalization_value = pure_normalized_value;
metadata.density_silent_normalization_value = density_normalized_value;
metadata.negative_prior_accepted_value = negative_prior_value;
metadata.unequal_pure_prior_formula_defect = ...
    ['For unequal priors, the pinned two-pure-state scalar uses ', ...
     'sqrt(2*(p1^2+p2^2)-4*p1*p2*overlap^2). The returned Helstrom POVM ', ...
     'attains the correct smaller value, so the pinned scalar can exceed 1.'];
metadata.evidence_scope = ...
    ['Deterministic solver-free one-state, two-state Helstrom, and exactly ', ...
     'orthogonal branches only. CVX was not available, so three-state ', ...
     'nonorthogonal SDP evidence is supplied independently by Hypatia/SCS. ', ...
     'The unequal-prior pure-state formula defect is recorded and corrected. ', ...
     'Silent state normalization and acceptance of a negative prior are ', ...
     'recorded as unsafe pinned behaviors and are deliberately rejected by ', ...
     'the Julia APIs.'];
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
fprintf('Wrote %d Distinguishability fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, input_kind, states, priors)
    [value, measurement] = Distinguishability(states, priors);
    if iscell(states)
        dimension = size(states{1}, 1);
        state_count = length(states);
        state_data = cat(3, states{:});
        normalized_states = cell(size(states));
        for index = 1:state_count
            normalized_states{index} = states{index} / trace(states{index});
        end
    else
        [dimension, state_count] = size(states);
        state_data = states;
        normalized_columns = normalize_cols(states);
        normalized_states = cell(1, state_count);
        for index = 1:state_count
            normalized_states{index} = ...
                normalized_columns(:, index) * normalized_columns(:, index)';
        end
    end
    if iscell(measurement)
        measurement_data = cat(3, measurement{:});
    else
        measurement_data = reshape(measurement, dimension, dimension, 1);
    end
    measurement_sum = sum(measurement_data, 3);
    completeness_residual = max(max(abs(measurement_sum - eye(dimension))));
    minimum_eigenvalue = Inf;
    objective = 0;
    for index = 1:state_count
        effect = (measurement_data(:, :, index) + ...
            measurement_data(:, :, index)') / 2;
        minimum_eigenvalue = min(minimum_eigenvalue, min(eig(effect)));
        objective = objective + priors(index) * ...
            real(trace(measurement_data(:, :, index) * normalized_states{index}));
    end

    fixture = struct();
    fixture.name = name;
    fixture.input_kind = input_kind;
    fixture.dimension = dimension;
    fixture.state_count = state_count;
    fixture.state_dims = size(state_data);
    fixture.state_real = reshape(full(real(state_data)), 1, []);
    fixture.state_imaginary = reshape(full(imag(state_data)), 1, []);
    fixture.priors = reshape(priors, 1, []);
    fixture.expected_value = value;
    fixture.measurement_real = reshape(full(real(measurement_data)), 1, []);
    fixture.measurement_imaginary = reshape(full(imag(measurement_data)), 1, []);
    fixture.measurement_completeness_residual = completeness_residual;
    fixture.measurement_positivity_violation = max(0, -minimum_eigenvalue);
    fixture.measurement_objective = objective;
    fixture.measurement_objective_residual = abs(objective - value);
    fixture.atol = 2e-10;
    fixture.rtol = 2e-10;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

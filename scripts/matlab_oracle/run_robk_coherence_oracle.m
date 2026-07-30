function run_robk_coherence_oracle()
% RUN_ROBK_COHERENCE_ORACLE Generate deterministic QETLAB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% RobkCohValue.m source SHA-256:
% 99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5
% QETLAB is Copyright 2014 Nathaniel Johnston and named coauthors,
% BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
if ~exist(fullfile(qetlab_path, 'RobkCohValue.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain RobkCohValue.m.');
end

addpath(qetlab_path);
cleanup_path = onCleanup(@() rmpath(qetlab_path));

fixtures = struct( ...
    'name', {}, ...
    'state_real', {}, ...
    'state_imaginary', {}, ...
    'k', {}, ...
    'robustness_real', {}, ...
    'robustness_imaginary', {}, ...
    'branch_index', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

fixtures(end + 1) = oracle_fixture( ...
    'basis_k2', [1; 0; 0; 0], 2);
fixtures(end + 1) = oracle_fixture( ...
    'rank_two_k2', normalized([1; 1; 0; 0]), 2);
fixtures(end + 1) = oracle_fixture( ...
    'uniform_four_k2', normalized([1; 1; 1; 1]), 2);
fixtures(end + 1) = oracle_fixture( ...
    'uniform_four_k3', normalized([1; 1; 1; 1]), 3);
fixtures(end + 1) = oracle_fixture( ...
    'strict_middle_branch_k3', normalized([0.7; 0.5; 0.4; sqrt(0.1)]), 3);
fixtures(end + 1) = oracle_fixture( ...
    'equality_branch_two_k3', normalized([0.5; 0.375; 0.375; 0.25]), 3);
fixtures(end + 1) = oracle_fixture( ...
    'equality_branch_three_k3', normalized([0.75; 0.5; 0.375; 0.125]), 3);
fixtures(end + 1) = oracle_fixture( ...
    'generic_five_k4', normalized([0.61; 0.47; 0.38; 0.31; 0.19]), 4);

canonical = normalized([0.7; 0.5; 0.4; sqrt(0.1)]);
correction_cases = struct( ...
    'name', {}, ...
    'state_real', {}, ...
    'state_imaginary', {}, ...
    'canonical_magnitudes', {}, ...
    'k', {}, ...
    'qetlab_robustness_real', {}, ...
    'qetlab_robustness_imaginary', {}, ...
    'qetlab_branch_index', {}, ...
    'expected_native_behavior', {});
correction_cases(end + 1) = correction_fixture( ...
    'unsorted_coefficients', canonical([4, 1, 3, 2]), canonical, 3, ...
    'sort_magnitudes_before_applying_the_theorem');
phased = canonical .* exp(1i * [0.3; -1.1; 2.4; 0.8]);
correction_cases(end + 1) = correction_fixture( ...
    'complex_phases', phased, canonical, 3, ...
    'discard_basis_phases_via_coefficient_magnitudes');
correction_cases(end + 1) = correction_fixture( ...
    'nonnormalized_coefficients', 2 * canonical, 2 * canonical, 3, ...
    'reject_without_implicit_normalization');

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-robk-coherence';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_source_sha256 = ...
    '99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5';
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.correction_case_count = length(correction_cases);
metadata.source_free_fixture = true;
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
payload.correction_cases = correction_cases;
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
fprintf('Wrote %d RobkCohValue fixtures and %d correction cases to %s using %s %s.\n', ...
    length(fixtures), length(correction_cases), output_path, ...
    metadata.engine, metadata.engine_version);
end

function vector = normalized(vector)
    vector = vector / norm(vector);
end

function fixture = oracle_fixture(name, state, k)
    [robustness, branch_index] = RobkCohValue(state, k);
    fixture = struct();
    fixture.name = name;
    fixture.state_real = reshape(full(real(state)), 1, []);
    fixture.state_imaginary = reshape(full(imag(state)), 1, []);
    fixture.k = k;
    fixture.robustness_real = real(robustness);
    fixture.robustness_imaginary = imag(robustness);
    fixture.branch_index = branch_index;
    fixture.comparison = 'approximate';
    fixture.atol = 5e-14;
    fixture.rtol = 5e-14;
end

function fixture = correction_fixture(name, state, canonical_magnitudes, k, behavior)
    [robustness, branch_index] = RobkCohValue(state, k);
    fixture = struct();
    fixture.name = name;
    fixture.state_real = reshape(full(real(state)), 1, []);
    fixture.state_imaginary = reshape(full(imag(state)), 1, []);
    fixture.canonical_magnitudes = ...
        reshape(full(real(canonical_magnitudes)), 1, []);
    fixture.k = k;
    fixture.qetlab_robustness_real = real(robustness);
    fixture.qetlab_robustness_imaginary = imag(robustness);
    fixture.qetlab_branch_index = branch_index;
    fixture.expected_native_behavior = behavior;
end

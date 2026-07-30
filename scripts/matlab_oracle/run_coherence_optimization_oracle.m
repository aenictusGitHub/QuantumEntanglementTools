function run_coherence_optimization_oracle()
% RUN_COHERENCE_OPTIMIZATION_ORACLE Generate deterministic analytic fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014 Nathaniel Johnston and named coauthors,
% BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.
%
% Every fixture below terminates in a solver-free theorem branch. The
% generalized-robustness fixtures use RobkCohValue, whose primary-source
% theorem gives the generalized and standard pure-state k-coherence
% robustness. No CVX result is recorded.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
required_files = {'IskIncoherent.m', 'IsAbskIncoh.m', ...
    'RobustnessCoherence.m', 'TraceDistanceCoherence.m', ...
    'GenRobustnesskCoherence.m', 'RobkCohValue.m'};
for index = 1:length(required_files)
    if ~exist(fullfile(qetlab_path, required_files{index}), 'file')
        error('QuantumEntanglementTools:OracleConfiguration', ...
            'The selected QETLAB checkout is missing %s.', required_files{index});
    end
end
helpers_path = fullfile(qetlab_path, 'helpers');
if ~exist(fullfile(helpers_path, 'has_band_k_ordering.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout is missing has_band_k_ordering.m.');
end

addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

isk_fixtures = struct('name', {}, 'state', {}, 'k', {}, 'qetlab_value', {});
isk_fixtures(end + 1) = isk_fixture( ...
    'diagonal_level_one', diag([0.4, 0.35, 0.25]), 1);
psi3 = ones(3, 1) / sqrt(3);
isk_fixtures(end + 1) = isk_fixture( ...
    'maximally_coherent_three_level_two', psi3 * psi3', 2);
isk_fixtures(end + 1) = isk_fixture( ...
    'full_dimension', psi3 * psi3', 3);
comparison_state = [0.40, 0.10, 0.00; 0.10, 0.35, 0.05; 0.00, 0.05, 0.25];
isk_fixtures(end + 1) = isk_fixture( ...
    'comparison_matrix_certificate', comparison_state, 2);
psi4 = ones(4, 1) / 2;
dephasing_state = 0.5 * eye(4) / 4 + 0.5 * (psi4 * psi4');
isk_fixtures(end + 1) = isk_fixture( ...
    'dephasing_certificate', dephasing_state, 3);

absolute_fixtures = struct('name', {}, 'state', {}, 'k', {}, 'qetlab_value', {});
absolute_fixtures(end + 1) = absolute_fixture( ...
    'maximally_mixed_level_one', eye(4) / 4, 1);
absolute_fixtures(end + 1) = absolute_fixture( ...
    'pure_rank_failure', diag([1, 0, 0, 0]), 2);
absolute_fixtures(end + 1) = absolute_fixture( ...
    'tight_rank_equal_spectrum', diag([1 / 3, 1 / 3, 1 / 3, 0]), 2);
absolute_fixtures(end + 1) = absolute_fixture( ...
    'low_dimension_purity_failure', diag([0.7, 0.2, 0.1]), 2);
absolute_fixtures(end + 1) = absolute_fixture( ...
    'one_sided_unknown', diag([0.55, 0.18, 0.12, 0.09, 0.06]), 3);

robustness_fixtures = struct( ...
    'name', {}, 'state', {}, 'qetlab_value', {}, 'atol', {}, 'rtol', {});
robustness_fixtures(end + 1) = robustness_fixture( ...
    'uniform_pure_four', psi4);
robustness_fixtures(end + 1) = robustness_fixture( ...
    'qubit_mixed', [0.6, 0.2i; -0.2i, 0.4]);
robustness_fixtures(end + 1) = robustness_fixture( ...
    'basis_state', [1; 0; 0]);

trace_fixtures = struct( ...
    'name', {}, 'state', {}, 'qetlab_value', {}, ...
    'closest_diagonal', {}, 'qetlab_output_dims', {}, 'atol', {}, 'rtol', {});
trace_fixtures(end + 1) = trace_fixture('uniform_pure_four', psi4);
trace_fixtures(end + 1) = trace_fixture( ...
    'qubit_mixed', [0.6, 0.2i; -0.2i, 0.4]);
generic_pure = normalized([0.7; 0.5; 0.4; sqrt(0.1)]);
trace_fixtures(end + 1) = trace_fixture('generic_pure_four', generic_pure);
trace_fixtures(end + 1) = trace_fixture('basis_state', [1; 0; 0]);

generalized_fixtures = struct( ...
    'name', {}, 'state', {}, 'k', {}, 'qetlab_value', {}, ...
    'qetlab_branch_index', {}, 'atol', {}, 'rtol', {});
generalized_fixtures(end + 1) = generalized_fixture( ...
    'uniform_four_level_two', psi4, 2);
generalized_fixtures(end + 1) = generalized_fixture( ...
    'uniform_four_level_three', psi4, 3);
generalized_fixtures(end + 1) = generalized_fixture( ...
    'generic_four_level_three', generic_pure, 3);

% Minimal counterexample found by exhaustive enumeration of all labelled
% three-vertex graphs. Its one-based bandwidth is 2 under ordering [2,1,3].
% Pinned QETLAB returns false because add_node calls order_is_reversed with
% candidate and num_placed swapped.
band_matrix = [0, 1, 1; 1, 0, 0; 1, 0, 0];
[pinned_has_ordering, pinned_ordering] = has_band_k_ordering(band_matrix, 2);
band_defect = struct();
band_defect.matrix = encode_array(band_matrix);
band_defect.k = 2;
band_defect.pinned_has_ordering = logical(pinned_has_ordering);
band_defect.pinned_ordering = reshape(pinned_ordering, 1, []);
band_defect.bruteforce_has_ordering = true;
band_defect.valid_ordering = [2, 1, 3];
band_defect.expected_native_disposition = ...
    'correct_argument_order_and_use_bounded_exact_layout_search';

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP4-coherence-optimization';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.source_free_fixture = true;
metadata.solver_free_branches_only = true;
metadata.isk_source_sha256 = ...
    'd924a12b239a1ea9a784ed3252ff0f289f2281aeae6022b645b169c4d2779f6a';
metadata.is_absk_source_sha256 = ...
    '6bf68ba8b87e441dd6b2ab2b718f621e2fa6a8a67a92f016b3629266496658ef';
metadata.robustness_source_sha256 = ...
    'b8af3c8f056ba49ad703e361b1b6f8554776747b0790efefed34a359a28b7986';
metadata.trace_distance_source_sha256 = ...
    'f74de030eb833da8c1bd2352961e1080c652746fa11572812000f50a0cd2048b';
metadata.generalized_source_sha256 = ...
    'aa3d429fc0fb9f715391620de67d9407f638bd58307922e49bf7ab5b57e37e78';
metadata.band_helper_source_sha256 = ...
    'd4b6b6be5793515b3b98f98ba585f04041e6e7875f22437a36b413f6308be8d7';
metadata.robk_source_sha256 = ...
    '99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5';
metadata.isk_fixture_count = length(isk_fixtures);
metadata.absolute_fixture_count = length(absolute_fixtures);
metadata.robustness_fixture_count = length(robustness_fixtures);
metadata.trace_fixture_count = length(trace_fixtures);
metadata.generalized_fixture_count = length(generalized_fixtures);
metadata.trace_qubit_output_defect = ...
    ['Pinned TraceDistanceCoherence documents a vector of diagonal entries ', ...
     'but returns a diagonal matrix on its qubit branch.'];
metadata.generalized_dependency_disposition = ...
    ['Pinned GenRobustnesskCoherence calls absent IskCoherent.m. Pure-state ', ...
     'fixtures therefore use the equal generalized/standard theorem in ', ...
     'pinned RobkCohValue; mixed-state SDP evidence is independent.'];
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
payload.isk_fixtures = isk_fixtures;
payload.absolute_fixtures = absolute_fixtures;
payload.robustness_fixtures = robustness_fixtures;
payload.trace_fixtures = trace_fixtures;
payload.generalized_fixtures = generalized_fixtures;
payload.band_defect = band_defect;
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
fprintf(['Wrote %d k-incoherence, %d absolute, %d robustness, %d trace, ', ...
    'and %d generalized fixtures to %s using %s %s.\n'], ...
    length(isk_fixtures), length(absolute_fixtures), ...
    length(robustness_fixtures), length(trace_fixtures), ...
    length(generalized_fixtures), output_path, ...
    metadata.engine, metadata.engine_version);
end

function fixture = isk_fixture(name, state, k)
fixture = struct();
fixture.name = name;
fixture.state = encode_array(state);
fixture.k = k;
fixture.qetlab_value = IskIncoherent(state, k);
end

function fixture = absolute_fixture(name, state, k)
fixture = struct();
fixture.name = name;
fixture.state = encode_array(state);
fixture.k = k;
fixture.qetlab_value = IsAbskIncoh(state, k);
end

function fixture = robustness_fixture(name, state)
fixture = struct();
fixture.name = name;
fixture.state = encode_array(state);
fixture.qetlab_value = real(RobustnessCoherence(state));
fixture.atol = 8e-13;
fixture.rtol = 8e-13;
end

function fixture = trace_fixture(name, state)
[value, diagonal_output] = TraceDistanceCoherence(state);
if isvector(diagonal_output)
    closest_diagonal = reshape(diagonal_output, 1, []);
else
    closest_diagonal = reshape(diag(diagonal_output), 1, []);
end
fixture = struct();
fixture.name = name;
fixture.state = encode_array(state);
fixture.qetlab_value = real(value);
fixture.closest_diagonal = closest_diagonal;
fixture.qetlab_output_dims = size(diagonal_output);
fixture.atol = 8e-13;
fixture.rtol = 8e-13;
end

function fixture = generalized_fixture(name, state, k)
[value, branch_index] = RobkCohValue(state, k);
fixture = struct();
fixture.name = name;
fixture.state = encode_array(state);
fixture.k = k;
fixture.qetlab_value = real(value);
fixture.qetlab_branch_index = branch_index;
fixture.atol = 8e-13;
fixture.rtol = 8e-13;
end

function encoded = encode_array(value)
encoded = struct();
encoded.dims = size(value);
encoded.real = reshape(full(real(value)), 1, []);
encoded.imaginary = reshape(full(imag(value)), 1, []);
end

function vector = normalized(vector)
vector = vector / norm(vector);
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

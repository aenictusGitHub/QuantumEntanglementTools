function run_nonlocal_games_oracle()
% RUN_NONLOCAL_GAMES_ORACLE Generate deterministic solver-free QETLAB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2014-2022 Nathaniel Johnston, Mateus Araujo,
% and Vincent Russo, BSD-2-Clause. Full upstream terms:
% licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
for source = {'XORGameValue.m', 'BellInequalityMax.m'}
    if ~exist(fullfile(qetlab_path, source{1}), 'file')
        error('QuantumEntanglementTools:OracleConfiguration', ...
            'The selected QETLAB checkout does not contain %s.', source{1});
    end
end
helpers_path = fullfile(qetlab_path, 'helpers');
if ~exist(fullfile(helpers_path, 'bcs_to_nonlocal.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain bcs_to_nonlocal.m.');
end

addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

xor_fixtures = struct( ...
    'name', {}, ...
    'probability_dims', {}, ...
    'probabilities', {}, ...
    'parity_dims', {}, ...
    'parity', {}, ...
    'expected', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

probabilities = ones(2, 2) / 4;
parity = [0, 0; 0, 1];
xor_fixtures(end + 1) = xor_fixture( ...
    'xor_chsh_classical', probabilities, parity, ...
    XORGameValue(probabilities, parity, 'classical'));

probabilities = [1, 2, 3; 4, 5, 6] / 21;
parity = [0, 0, 0; 0, 1, 1];
xor_fixtures(end + 1) = xor_fixture( ...
    'xor_rectangular_classical', probabilities, parity, ...
    XORGameValue(probabilities, parity, 'classical'));

bell_fixtures = struct( ...
    'name', {}, ...
    'desc', {}, ...
    'notation', {}, ...
    'coefficient_dims', {}, ...
    'coefficients', {}, ...
    'expected', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

coefficients = [0, 0, 0; 0, 1, 1; 0, 1, -1];
bell_fixtures(end + 1) = bell_fixture( ...
    'bell_chsh_full_correlator_classical', ...
    coefficients, [2, 2, 2, 2], 'fc', ...
    BellInequalityMax(coefficients, [2, 2, 2, 2], 'fc', 'classical'));

coefficients = [3, 2, -1; -2, 1, 1; 1, 1, -1];
bell_fixtures(end + 1) = bell_fixture( ...
    'bell_affine_full_correlator_classical', ...
    coefficients, [2, 2, 2, 2], 'fc', ...
    BellInequalityMax(coefficients, [2, 2, 2, 2], 'fc', 'classical'));

coefficients = reshape(-5:18, [3, 2, 2, 2]);
bell_fixtures(end + 1) = bell_fixture( ...
    'bell_ternary_binary_full_probability_classical', ...
    coefficients, [3, 2, 2, 2], 'fp', ...
    BellInequalityMax(coefficients, [3, 2, 2, 2], 'fp', 'classical'));

constraints = {[0, 1; 1, 0], [1, 0; 1, 0]};
[bcs_probabilities, bcs_payoff] = bcs_to_nonlocal(constraints);
bcs_fixture = struct();
bcs_fixture.name = 'bcs_to_nonlocal_active_and_inactive_variables';
bcs_fixture.constraint_1_dims = size(constraints{1});
bcs_fixture.constraint_1 = reshape(constraints{1}, 1, []);
bcs_fixture.constraint_2_dims = size(constraints{2});
bcs_fixture.constraint_2 = reshape(constraints{2}, 1, []);
bcs_fixture.probability_dims = size(bcs_probabilities);
bcs_fixture.probabilities = reshape(bcs_probabilities, 1, []);
bcs_fixture.payoff_dims = size(bcs_payoff);
bcs_fixture.payoff = reshape(bcs_payoff, 1, []);
bcs_fixture.comparison = 'exact';
bcs_fixture.atol = 0;
bcs_fixture.rtol = 0;

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP6-nonlocal-games';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = ...
    length(xor_fixtures) + length(bell_fixtures) + 1;
metadata.source_free_fixture = true;
metadata.executed_scope = [ ...
    'solver-free classical XORGameValue and BellInequalityMax branches; ', ...
    'bcs_to_nonlocal helper'];
metadata.excluded_scope = [ ...
    'CVX-dependent NPA, no-signalling, quantum, fixed-qubit, and lower-bound ', ...
    'branches were not executed'];
metadata.upstream_blocker = [ ...
    'BCSGameValue calls NonlocalGameValue.m, which is absent from the pinned ', ...
    'QETLAB checkout'];
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
payload.xor_fixtures = xor_fixtures;
payload.bell_fixtures = bell_fixtures;
payload.bcs_fixture = bcs_fixture;
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
fprintf('Wrote %d nonlocal-game fixtures to %s using %s %s.\n', ...
    metadata.fixture_count, output_path, metadata.engine, metadata.engine_version);
end

function fixture = xor_fixture(name, probabilities, parity, expected)
    fixture = struct();
    fixture.name = name;
    fixture.probability_dims = size(probabilities);
    fixture.probabilities = reshape(probabilities, 1, []);
    fixture.parity_dims = size(parity);
    fixture.parity = reshape(parity, 1, []);
    fixture.expected = expected;
    fixture.comparison = 'approximate';
    fixture.atol = 1e-14;
    fixture.rtol = 1e-14;
end

function fixture = bell_fixture(name, coefficients, desc, notation, expected)
    fixture = struct();
    fixture.name = name;
    fixture.desc = desc;
    fixture.notation = notation;
    fixture.coefficient_dims = size(coefficients);
    fixture.coefficients = reshape(coefficients, 1, []);
    fixture.expected = expected;
    fixture.comparison = 'exact';
    fixture.atol = 0;
    fixture.rtol = 0;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

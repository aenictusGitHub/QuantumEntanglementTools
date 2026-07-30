function run_random_superoperator_oracle()
% RUN_RANDOM_SUPEROPERATOR_ORACLE Generate seeded QETLAB property fixtures.
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
if ~exist(fullfile(qetlab_path, 'RandomSuperoperator.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain RandomSuperoperator.m.');
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

fixtures = struct( ...
    'name', {}, ...
    'seed', {}, ...
    'dimensions', {}, ...
    'trace_preserving', {}, ...
    'unital', {}, ...
    'real_output', {}, ...
    'kraus_rank', {}, ...
    'matrix_dimensions', {}, ...
    'real', {}, ...
    'imaginary', {}, ...
    'numerical_rank', {}, ...
    'minimum_eigenvalue', {}, ...
    'choi_trace', {}, ...
    'trace_preservation_residual', {}, ...
    'unitality_residual', {}, ...
    'proportional_unitality_residual', {});

fixtures(end + 1) = seeded_fixture( ...
    'unconstrained_complex', 4101, [2, 3], 0, 0, 0, 3);
fixtures(end + 1) = seeded_fixture( ...
    'trace_preserving_real', 4102, [3, 3], 1, 0, 1, 2);
fixtures(end + 1) = seeded_fixture( ...
    'unital_real', 4103, [3, 3], 0, 1, 1, 2);
fixtures(end + 1) = seeded_fixture( ...
    'bistochastic_rank_one_real', 4104, [2, 2], 1, 1, 1, 1);

unequal_balance_error_message = '';
if exist('OCTAVE_VERSION', 'builtin')
    try
        fixtures(end + 1) = seeded_fixture( ...
            'unequal_balanced_is_not_unital', 4105, [2, 3], 1, 1, 1, 2);
    catch unequal_error
        unequal_balance_error_message = unequal_error.message;
    end
else
    fixtures(end + 1) = seeded_fixture( ...
        'unequal_balanced_is_not_unital', 4105, [2, 3], 1, 1, 1, 2);
end

rng(4199, 'twister');
permissive_flag_output = RandomSuperoperator([2, 3], 2, 0, 0, 2);
permissive_flag_trace = trace(permissive_flag_output);
rng(4200, 'twister');
oversized_rank_output = RandomSuperoperator([2, 2], 0, 0, 0, 5);
oversized_rank_numerical_rank = rank(oversized_rank_output, 1e-10);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-random-superoperator';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.source_sha256 = ...
    'ac2862092569e0878c0f5dac68f4cd58670992ae44adea26d2d7dfb36c07d96f';
metadata.comparison_scope = ...
    ['Seeded cross-engine draws are not compared entrywise. The artifact ', ...
     'records every executable TP/UN/RE branch and invariant residuals.'];
metadata.unequal_dimension_defect = ...
    ['TP=1 and UN=1 with unequal dimensions cannot be unital. The pinned ', ...
     'output maps identity to DIM(1)/DIM(2) times identity.'];
metadata.unequal_balance_error_message = unequal_balance_error_message;
metadata.octave_unequal_dimension_limitation = ...
    ['Octave 11.3 can fail inside the pinned PartialTrace cellfun path for ', ...
     'complex block data. Real constrained fixtures avoid treating that ', ...
     'engine-specific failure as MATLAB behavior.'];
metadata.permissive_flag_behavior = ...
    'The pinned routine accepts TP=2 and silently takes the unconstrained branch.';
metadata.permissive_flag_trace = permissive_flag_trace;
metadata.oversized_rank_behavior = ...
    ['The pinned routine accepts KR larger than prod(DIM), although the ', ...
     'documented almost-sure exact-rank statement is then impossible.'];
metadata.oversized_rank_numerical_rank = oversized_rank_numerical_rank;
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
fprintf('Wrote %d RandomSuperoperator property fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = seeded_fixture(name, seed, dimensions, tp, un, re, kr)
rng(seed, 'twister');
value = RandomSuperoperator(dimensions, tp, un, re, kr);
input_dimension = dimensions(1);
output_dimension = dimensions(2);
[input_marginal, output_marginal] = choi_marginals( ...
    value, input_dimension, output_dimension);
proportional_factor = input_dimension / output_dimension;

fixture = struct();
fixture.name = name;
fixture.seed = seed;
fixture.dimensions = dimensions;
fixture.trace_preserving = tp;
fixture.unital = un;
fixture.real_output = re;
fixture.kraus_rank = kr;
fixture.matrix_dimensions = size(value);
fixture.real = reshape(full(real(value)), 1, []);
fixture.imaginary = reshape(full(imag(value)), 1, []);
fixture.numerical_rank = rank(value, 1e-10);
fixture.minimum_eigenvalue = min(real(eig((value + value') / 2)));
fixture.choi_trace = real(trace(value));
fixture.trace_preservation_residual = ...
    norm(input_marginal - eye(input_dimension), 'fro');
fixture.unitality_residual = ...
    norm(output_marginal - eye(output_dimension), 'fro');
fixture.proportional_unitality_residual = ...
    norm(output_marginal - proportional_factor * eye(output_dimension), 'fro');
end

function [input_marginal, output_marginal] = choi_marginals( ...
    value, input_dimension, output_dimension)
input_marginal = zeros(input_dimension, input_dimension);
output_marginal = zeros(output_dimension, output_dimension);
for input_row = 1:input_dimension
    for input_column = 1:input_dimension
        for output_index = 1:output_dimension
            row = output_index + (input_row - 1) * output_dimension;
            column = output_index + (input_column - 1) * output_dimension;
            input_marginal(input_row, input_column) = ...
                input_marginal(input_row, input_column) + value(row, column);
        end
    end
end
for output_row = 1:output_dimension
    for output_column = 1:output_dimension
        for input_index = 1:input_dimension
            row = output_row + (input_index - 1) * output_dimension;
            column = output_column + (input_index - 1) * output_dimension;
            output_marginal(output_row, output_column) = ...
                output_marginal(output_row, output_column) + value(row, column);
        end
    end
end
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

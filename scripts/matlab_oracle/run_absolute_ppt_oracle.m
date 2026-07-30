function run_absolute_ppt_oracle()
% RUN_ABSOLUTE_PPT_ORACLE Generate deterministic QETLAB absolute-PPT fixtures.
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
required_files = {'AbsPPTConstraints.m', 'IsAbsPPT.m', 'IsPSD.m', ...
    'InSeparableBall.m'};
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

fixtures = struct( ...
    'name', {}, ...
    'dims', {}, ...
    'real', {}, ...
    'imaginary', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

lam2 = [0.4, 0.3, 0.2, 0.1];
constraints2 = AbsPPTConstraints(lam2, [2, 2]);
fixtures(end + 1) = oracle_fixture('p2_constraint_1', constraints2{1});

lam3 = [0.30, 0.19, 0.14, 0.11, 0.09, 0.07, 0.05, 0.03, 0.02];
constraints3 = AbsPPTConstraints(lam3, [3, 3]);
for index = 1:length(constraints3)
    fixtures(end + 1) = oracle_fixture( ...
        sprintf('p3_constraint_%d', index), constraints3{index});
end

lam4 = (16:-1:1) / sum(1:16);
constraints4 = AbsPPTConstraints(lam4, [4, 4], 0, 3);
for index = 1:length(constraints4)
    fixtures(end + 1) = oracle_fixture( ...
        sprintf('p4_limited_constraint_%d', index), constraints4{index});
end

pure4 = [1, zeros(1, 15)];
early4 = AbsPPTConstraints(pure4, [4, 4], 1);
fixtures(end + 1) = oracle_fixture('p4_early_violation', early4{end});

is_abs_ppt_outputs = struct();
is_abs_ppt_outputs.maximally_mixed_2x2 = IsAbsPPT(ones(1, 4) / 4, [2, 2]);
is_abs_ppt_outputs.exhaustive_yes_2x2 = ...
    IsAbsPPT([0.45, 0.35, 0.1, 0.1], [2, 2]);
is_abs_ppt_outputs.pure_no_2x2 = IsAbsPPT([1, 0, 0, 0], [2, 2]);
is_abs_ppt_outputs.boundary_2x2 = ...
    IsAbsPPT([0.5, 1 / 6, 1 / 6, 1 / 6], [2, 2]);
is_abs_ppt_outputs.default_rectangular_length6 = ...
    IsAbsPPT(ones(1, 6) / 6);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP3-absolute-ppt';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.abs_constraints_source_sha256 = ...
    'e7cdc38ce5b16a68165f0b983a6033496147fd8a6c1df4a33743c15bd11a5b75';
metadata.is_abs_ppt_source_sha256 = ...
    'f6272ff80fa53789ea92485f32b9155bf55327670d8e24c58f51f59f1286a5d2';
metadata.p2_constraint_count = length(constraints2);
metadata.p3_constraint_count = length(constraints3);
metadata.p4_full_constraint_count = length(AbsPPTConstraints(lam4, [4, 4]));
metadata.p4_limited_constraint_count = length(constraints4);
metadata.p4_early_constraint_count = length(early4);
metadata.scalar_dim_inconsistency = ...
    ['Pinned AbsPPTConstraints interprets scalar DIM as equal local dimensions, ', ...
     'while pinned IsAbsPPT interprets it as the first dimension and infers ', ...
     'the second. The native API uses the latter convention consistently.'];
metadata.p6_redundant_constraint_note = ...
    ['The pinned criss-cross enumeration returns 2612 constraints for p=6; ', ...
     'QETLAB documentation records four as redundant.'];
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
payload.is_abs_ppt_outputs = is_abs_ppt_outputs;
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
fprintf('Wrote %d absolute-PPT fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value)
fixture = struct();
fixture.name = name;
fixture.dims = size(value);
fixture.real = reshape(full(real(value)), 1, []);
fixture.imaginary = reshape(full(imag(value)), 1, []);
fixture.comparison = 'normwise';
fixture.atol = 8e-14;
fixture.rtol = 8e-14;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
rmpath(helpers_path);
rmpath(qetlab_path);
end

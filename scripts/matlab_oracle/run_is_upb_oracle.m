function run_is_upb_oracle()
% RUN_IS_UPB_ORACLE Generate deterministic pinned-QETLAB IsUPB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB IsUPB.m and helpers/vec_partitions.m are Copyright 2014
% Nathaniel Johnston, BSD-2-Clause.
% Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
if ~exist(fullfile(qetlab_path, 'IsUPB.m'), 'file') || ...
        ~exist(fullfile(qetlab_path, 'helpers', 'vec_partitions.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout lacks IsUPB.m or vec_partitions.m.');
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

tiles_left = [ ...
    1, 1, 0, 0, 1; ...
    0, -1, 0, 1, 1; ...
    0, 0, 1, -1, 1];
tiles_right = [ ...
    1, 0, 0, 1, 1; ...
    -1, 0, 1, 0, 1; ...
    0, 1, -1, 0, 1];
tiles_flag = IsUPB(tiles_left, tiles_right);

shifts_first = [1, 0, 1, 1; 0, 1, -1, 1];
shifts_second = [1, 1, 0, 1; 0, 1, 1, -1];
shifts_third = [1, 1, 1, 0; 0, -1, 1, 1];
shifts_flag = IsUPB(shifts_first, shifts_second, shifts_third);

extendible_left = [1, 0; 0, 1];
extendible_right = [1, 1; 0, 0];
[extendible_flag, extendible_witness] = ...
    IsUPB(extendible_left, extendible_right);

complete_left = [1, 1, 0, 0; 0, 0, 1, 1];
complete_right = [1, 0, 1, 0; 0, 1, 0, 1];
complete_flag = IsUPB(complete_left, complete_right);

generic_left = [ ...
    1, 1, 1, 1, 1; ...
    0, 1, 2, 3, 4; ...
    0, 1, 4, 9, 16];
generic_right = [ ...
    1, 1, 1, 1, 1; ...
    0, 2, 4, 6, 8; ...
    0, 4, 16, 36, 64];
nonorthogonal_flag = IsUPB(generic_left, generic_right);
nonorthogonal_overlap = ...
    abs((generic_left(:, 1)' * generic_left(:, 2)) * ...
        (generic_right(:, 1)' * generic_right(:, 2)));

complex_left = [1; 0];
complex_right = [1; 1i];
[complex_flag, complex_witness] = IsUPB(complex_left, complex_right);
complex_hilbert_residual = ...
    abs((complex_witness{1}' * complex_left) * ...
        (complex_witness{2}' * complex_right));
complex_bilinear_residual = ...
    abs((complex_witness{1}.' * complex_left) * ...
        (complex_witness{2}.' * complex_right));
conjugated_hilbert_residual = ...
    abs((conj(complex_witness{1})' * complex_left) * ...
        (conj(complex_witness{2})' * complex_right));

fixtures(end + 1) = oracle_fixture('is_upb_tiles_flag', tiles_flag);
fixtures(end + 1) = oracle_fixture('is_upb_shifts_flag', shifts_flag);
fixtures(end + 1) = oracle_fixture('is_upb_extendible_flag', extendible_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_extendible_witness_left', extendible_witness{1});
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_extendible_witness_right', extendible_witness{2});
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_complete_basis_pinned_flag', complete_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_nonorthogonal_pinned_flag', nonorthogonal_flag);
fixtures(end + 1) = oracle_fixture('is_upb_complex_flag', complex_flag);
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_complex_witness_left', complex_witness{1});
fixtures(end + 1) = oracle_fixture( ...
    'is_upb_complex_witness_right', complex_witness{2});

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-is-upb';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.nonorthogonal_first_overlap = nonorthogonal_overlap;
metadata.complex_witness_hilbert_residual = complex_hilbert_residual;
metadata.complex_witness_bilinear_residual = complex_bilinear_residual;
metadata.conjugated_witness_hilbert_residual = conjugated_hilbert_residual;
metadata.reviewed_deviations = [ ...
    'The pinned routine checks extendibility but not mutual orthogonality ', ...
    'or incompleteness, and its nonconjugating-transpose witness need not ', ...
    'be Hilbert-space orthogonal for complex local factors.'];
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
fprintf('Wrote %d IsUPB fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = oracle_fixture(name, value)
    fixture = struct();
    fixture.name = name;
    fixture.dims = size(value);
    fixture.real = reshape(full(real(value)), 1, []);
    fixture.imaginary = reshape(full(imag(value)), 1, []);
    fixture.comparison = 'exact';
    fixture.atol = 0;
    fixture.rtol = 0;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

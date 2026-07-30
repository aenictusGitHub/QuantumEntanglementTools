function run_upb_catalog_oracle()
% RUN_UPB_CATALOG_ORACLE Generate deterministic pinned-QETLAB UPB fixtures.
%
% QETLAB source revision:
% d8589610f00cff106537268dee2e2a1153f3a601
% QETLAB is Copyright 2022 Nathaniel Johnston, BSD-2-Clause.
% Full upstream terms: licenses/QETLAB-LICENSE.txt.

qetlab_path = getenv('QET_ORACLE_QETLAB_PATH');
output_path = getenv('QET_ORACLE_OUTPUT');
qetlab_commit = getenv('QET_ORACLE_QETLAB_COMMIT');
if isempty(qetlab_path) || isempty(output_path) || isempty(qetlab_commit)
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'QETLAB path, output path, and commit environment variables are required.');
end
if ~exist(fullfile(qetlab_path, 'UPB.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain UPB.m.');
end

helpers_path = fullfile(qetlab_path, 'helpers');
addpath(qetlab_path);
addpath(helpers_path);
cleanup_path = onCleanup(@() remove_qetlab_paths(qetlab_path, helpers_path));

fixtures = struct( ...
    'name', {}, ...
    'family', {}, ...
    'dimensions', {}, ...
    'cardinality', {}, ...
    'factor_rows', {}, ...
    'real', {}, ...
    'imaginary', {}, ...
    'comparison', {}, ...
    'atol', {}, ...
    'rtol', {});

fixtures(end + 1) = named_fixture('pyramid', 'Pyramid', {}, 2);
fixtures(end + 1) = named_fixture('tiles', 'Tiles', {}, 2);
fixtures(end + 1) = named_fixture('generalized_tiles_1_d4', ...
    'GenTiles1', {4}, 2);
fixtures(end + 1) = named_fixture('generalized_tiles_2_3x4', ...
    'GenTiles2', {3, 4}, 2);
fixtures(end + 1) = named_fixture('minimum_4x4', 'Min4x4', {}, 2);
fixtures(end + 1) = named_fixture('quad_residue_d3', 'QuadRes', {3}, 2);
fixtures(end + 1) = named_fixture('six_parameter', 'SixParam', ...
    {[0.31, 0.47, -0.2, 0.53, 0.61, 0.4]}, 2);
fixtures(end + 1) = named_fixture('generalized_shifts_p5', ...
    'GenShifts', {5}, 5);
fixtures(end + 1) = named_fixture('feng_2x2x2x2', ...
    'Feng2x2x2x2', {}, 4);
fixtures(end + 1) = named_fixture('johnston_2_power_8', ...
    'John2^8', {}, 8);
fixtures(end + 1) = named_fixture('chen_johnston_4x6', ...
    'CJBip46', {}, 2);
fixtures(end + 1) = named_fixture('feng_2x2x3', ...
    'Feng2x2x3', {}, 3);
fixtures(end + 1) = named_fixture('feng_2x2x5', ...
    'Feng2x2x5', {}, 3);
fixtures(end + 1) = named_fixture('feng_4x4', 'Feng4x4', {}, 2);
fixtures(end + 1) = named_fixture('feng_2x2x2x4', ...
    'Feng2x2x2x4', {}, 4);
fixtures(end + 1) = named_fixture('feng_2x2x2x2x5', ...
    'Feng2x2x2x2x5', {}, 5);
fixtures(end + 1) = dimension_fixture('dimension_dispatch_3x2', [3, 2]);

defect_locals = cell(1, 8);
[defect_locals{:}] = UPB('John2^4k', 8);
defect_gram = ones(12);
for party = 1:length(defect_locals)
    defect_gram = defect_gram .* (defect_locals{party}' * defect_locals{party});
end
defect = struct();
defect.name = 'johnston_2_power_4k_reshape_order';
defect.family = 'John2^4k';
defect.parameter = 8;
defect.cardinality = 12;
defect.pinned_orthogonality_residual = ...
    max(max(abs(defect_gram - eye(12))));
defect.native_disposition = ...
    ['corrected from the primary-source orthogonality graph; ', ...
     'entrywise parity is intentionally not required'];

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'WP2-UPB-catalog';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.upb_sha256 = ...
    'd60d9d2375d773c7e59930991ddf73427c9ffb0c77ded640b7f635f99ec13544';
metadata.one_factorization_sha256 = ...
    '7b042b3729b96dc4646c57adc24dcd5e4f59cf710520b0a6f98e7f3fd0c531ca';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.source_free_fixture = true;
metadata.comparison_scope = ...
    ['deterministic executable families compared up to independent ', ...
     'local-vector phases; randomized families use theorem/property tests'];
metadata.upstream_defect_count = 1;
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
payload.upstream_defects = defect;
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
fprintf('Wrote %d UPB fixtures to %s using %s %s.\n', ...
    length(fixtures), output_path, metadata.engine, metadata.engine_version);
end

function fixture = named_fixture(name, family, arguments, party_count)
    local_factors = cell(1, party_count);
    [local_factors{:}] = UPB(family, arguments{:});
    fixture = local_fixture(name, family, local_factors);
end

function fixture = dimension_fixture(name, dimensions)
    local_factors = cell(1, length(dimensions));
    [local_factors{:}] = UPB(dimensions, 0);
    fixture = local_fixture(name, 'DIM', local_factors);
end

function fixture = local_fixture(name, family, local_factors)
    cardinality = size(local_factors{1}, 2);
    dimensions = zeros(1, length(local_factors));
    for party = 1:length(local_factors)
        dimensions(party) = size(local_factors{party}, 1);
        if size(local_factors{party}, 2) ~= cardinality
            error('QuantumEntanglementTools:OracleGeneration', ...
                'UPB local factors have inconsistent cardinalities.');
        end
    end
    concatenated = vertcat(local_factors{:});
    fixture = struct();
    fixture.name = name;
    fixture.family = family;
    fixture.dimensions = dimensions;
    fixture.cardinality = cardinality;
    fixture.factor_rows = size(concatenated, 1);
    fixture.real = reshape(real(concatenated), 1, []);
    fixture.imaginary = reshape(imag(concatenated), 1, []);
    fixture.comparison = 'local_vectors_up_to_independent_phase';
    fixture.atol = 2e-11;
    fixture.rtol = 2e-11;
end

function remove_qetlab_paths(qetlab_path, helpers_path)
    rmpath(helpers_path);
    rmpath(qetlab_path);
end

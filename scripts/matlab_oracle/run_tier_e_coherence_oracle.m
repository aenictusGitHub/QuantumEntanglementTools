function run_tier_e_coherence_oracle()
% RUN_TIER_E_COHERENCE_ORACLE Generate deterministic QETLAB coherence fixtures.
%
% The CoherenceRank fixtures deliberately record the pinned implementation's
% zero-counting bug. They are evidence of a reviewed upstream discrepancy, not
% expected values for the Julia-native mathematical coherence rank.
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
if ~exist(fullfile(qetlab_path, 'L1NormCoherence.m'), 'file')
    error('QuantumEntanglementTools:OracleConfiguration', ...
        'The selected QETLAB checkout does not contain L1NormCoherence.m.');
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

plus_ququart = ones(4, 1) / 2;
mixed_qubit = [0.6, 0.2i; -0.2i, 0.4];
plus_qubit = ones(2, 1) / sqrt(2);
entropy_state = [0.7, 0.1i; -0.1i, 0.3];

fixtures(end + 1) = oracle_fixture( ...
    'l1_plus_ququart', L1NormCoherence(plus_ququart), ...
    'normwise', 5e-13, 5e-13);
fixtures(end + 1) = oracle_fixture( ...
    'l1_mixed_qubit', L1NormCoherence(mixed_qubit), ...
    'normwise', 5e-13, 5e-13);
fixtures(end + 1) = oracle_fixture( ...
    'relative_entropy_plus_qubit_base2', RelEntCoherence(plus_qubit), ...
    'normwise', 5e-12, 5e-12);
fixtures(end + 1) = oracle_fixture( ...
    'relative_entropy_mixed_qubit_base2', RelEntCoherence(entropy_state), ...
    'normwise', 5e-12, 5e-12);

% QETLAB's implementation increments its counter for entries at or below the
% tolerance, contrary to its documentation's nonzero-entry definition.
fixtures(end + 1) = oracle_fixture( ...
    'coherence_rank_qetlab_plus_bug', CoherenceRank(plus_qubit), ...
    'exact', 0, 0);
fixtures(end + 1) = oracle_fixture( ...
    'coherence_rank_qetlab_basis_bug', CoherenceRank([1; 0; 0]), ...
    'exact', 0, 0);

metadata = struct();
metadata.schema = 'quantum-entanglement-tools-oracle-v1';
metadata.tier = 'E-coherence';
metadata.generated_utc = [datestr(now, 30), 'Z'];
metadata.qetlab_commit = qetlab_commit;
metadata.qetlab_license = 'BSD-2-Clause';
metadata.platform = computer;
metadata.fixture_count = length(fixtures);
metadata.reviewed_upstream_discrepancy = ...
    'CoherenceRank counts entries at or below tolerance instead of nonzero entries';

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

fprintf('Wrote %d Tier E coherence fixtures to %s using %s %s.\n', ...
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

#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "${script_directory}/../.." && pwd)"
qetlab_path="${repository_root}/dev/upstream/QETLAB"
output_path="${repository_root}/test/oracle/generated/symmetric_extensions_oracle.json"
engine="auto"
expected_commit="d8589610f00cff106537268dee2e2a1153f3a601"
expected_symmetric_extension_sha="ca8e1aaf9b3766a3e29bc1375eecf473cdb77bc71159a033411c77eacb6e7e59"
expected_symmetric_inner_extension_sha="f73e38499197fb18c1ff9cff554fd77e13ab7ccfbe48590a25437cd70d7f39da"
expected_random_ppt_state_sha="e5208eac2c04dac95051d8c5850bb6fdaaca6430117416d0675b2c9a26ef04fd"
expected_jacobi_poly_sha="3001148f6fb136306a23a48ce6176b2b2ed65b278e4f1f998cb180486124d513"

usage() {
    printf '%s\n' \
        "Usage: $0 [--engine auto|matlab|octave] [--qetlab PATH] [--output PATH]" \
        "" \
        "MATLAB is preferred in auto mode. Octave evidence is supplemental."
}

while (($# > 0)); do
    case "$1" in
        --engine)
            engine="${2:?--engine requires a value}"
            shift 2
            ;;
        --qetlab)
            qetlab_path="${2:?--qetlab requires a path}"
            shift 2
            ;;
        --output)
            output_path="${2:?--output requires a path}"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$engine" in
    auto|matlab|octave) ;;
    *)
        printf 'Invalid engine %s; expected auto, matlab, or octave.\n' "$engine" >&2
        exit 2
        ;;
esac

if [[ ! -f "${qetlab_path}/SymmetricExtension.m" ]]; then
    printf 'QETLAB checkout not found at %s.\n' "$qetlab_path" >&2
    exit 2
fi
actual_commit="$(git -C "$qetlab_path" rev-parse HEAD)"
if [[ "$actual_commit" != "$expected_commit" ]]; then
    printf 'Refusing unpinned QETLAB checkout: expected %s, found %s.\n' \
        "$expected_commit" "$actual_commit" >&2
    exit 2
fi

check_source_sha() {
    local relative_path="$1"
    local expected_sha="$2"
    local actual_sha
    actual_sha="$(shasum -a 256 "${qetlab_path}/${relative_path}" | awk '{print $1}')"
    if [[ "$actual_sha" != "$expected_sha" ]]; then
        printf 'Refusing %s with unexpected SHA-256: expected %s, found %s.\n' \
            "$relative_path" "$expected_sha" "$actual_sha" >&2
        exit 2
    fi
}

check_source_sha "SymmetricExtension.m" "$expected_symmetric_extension_sha"
check_source_sha "SymmetricInnerExtension.m" "$expected_symmetric_inner_extension_sha"
check_source_sha "RandomPPTState.m" "$expected_random_ppt_state_sha"
check_source_sha "helpers/jacobi_poly.m" "$expected_jacobi_poly_sha"

mkdir -p "$(dirname "$output_path")"
export QET_ORACLE_QETLAB_PATH="$qetlab_path"
export QET_ORACLE_OUTPUT="$output_path"
export QET_ORACLE_QETLAB_COMMIT="$actual_commit"

matlab_binary="${QET_ORACLE_MATLAB_BIN:-matlab}"
octave_binary="${QET_ORACLE_OCTAVE_BIN:-}"
if [[ -z "$octave_binary" ]]; then
    if command -v octave-cli >/dev/null 2>&1; then
        octave_binary="octave-cli"
    else
        octave_binary="octave"
    fi
fi
selected_engine="$engine"
if [[ "$selected_engine" == "auto" ]]; then
    if command -v "$matlab_binary" >/dev/null 2>&1; then
        selected_engine="matlab"
    elif command -v "$octave_binary" >/dev/null 2>&1; then
        selected_engine="octave"
    else
        printf 'Neither MATLAB nor Octave was found; oracle generation is optional.\n' >&2
        exit 2
    fi
fi

if [[ "$selected_engine" == "matlab" ]]; then
    command -v "$matlab_binary" >/dev/null 2>&1 || {
        printf 'MATLAB executable not found: %s\n' "$matlab_binary" >&2
        exit 2
    }
    escaped_directory="${script_directory//\'/\'\'}"
    "$matlab_binary" -batch \
        "addpath('${escaped_directory}'); run_symmetric_extensions_oracle"
else
    command -v "$octave_binary" >/dev/null 2>&1 || {
        printf 'Octave executable not found: %s\n' "$octave_binary" >&2
        exit 2
    }
    "$octave_binary" --quiet --no-gui --path "$script_directory" \
        --eval "run_symmetric_extensions_oracle"
fi

if command -v shasum >/dev/null 2>&1; then
    (
        cd "$(dirname "$output_path")"
        shasum -a 256 "$(basename "$output_path")" \
            > "$(basename "$output_path").sha256"
    )
elif command -v sha256sum >/dev/null 2>&1; then
    (
        cd "$(dirname "$output_path")"
        sha256sum "$(basename "$output_path")" \
            > "$(basename "$output_path").sha256"
    )
else
    printf 'No SHA-256 command found; fixture exists but hash was not recorded.\n' >&2
    exit 2
fi
printf 'Oracle fixture: %s\n' "$output_path"
printf 'Fixture digest: %s.sha256\n' "$output_path"

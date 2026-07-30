#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "${script_directory}/../.." && pwd)"
qetlab_path="${repository_root}/dev/upstream/QETLAB"
output_path="${repository_root}/test/oracle/generated/state_discrimination_oracle.json"
engine="auto"
expected_commit="d8589610f00cff106537268dee2e2a1153f3a601"

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

if [[ ! -f "${qetlab_path}/Distinguishability.m" ]]; then
    printf 'QETLAB checkout not found at %s.\n' "$qetlab_path" >&2
    exit 2
fi
actual_commit="$(git -C "$qetlab_path" rev-parse HEAD)"
if [[ "$actual_commit" != "$expected_commit" ]]; then
    printf 'Refusing unpinned QETLAB checkout: expected %s, found %s.\n' \
        "$expected_commit" "$actual_commit" >&2
    exit 2
fi

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
        "addpath('${escaped_directory}'); run_state_discrimination_oracle"
else
    command -v "$octave_binary" >/dev/null 2>&1 || {
        printf 'Octave executable not found: %s\n' "$octave_binary" >&2
        exit 2
    }
    "$octave_binary" --quiet --no-gui --path "$script_directory" \
        --eval "run_state_discrimination_oracle"
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


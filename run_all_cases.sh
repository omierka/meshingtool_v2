#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_all_cases.sh [-c|--clean]

Loops over CASES/* and runs runner.sh for each case. Use -c/--clean to add the
clean flag when invoking runner.sh. Executes runner.sh in silent mode (stage
markers only) and prints "case-name, timeconsumption [s]" for each case.
Use -n/--num-proc to forward a custom MPI process count to runner.sh (default 64).
EOF
}

num_proc="4"
clean_flag=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--clean)
      clean_flag+=("-c")
      shift
      ;;
    -n|--num-proc)
      shift
      if [[ $# -eq 0 ]]; then
        echo "Missing value for $1" >&2
        usage >&2
        exit 1
      fi
      num_proc="$1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cases_dir="${script_dir}/CASES"
runner="${script_dir}/newgenmeshingtool"
required_binaries=(
  "newgenmeshingtool"
  "meshhexer-cli"
  "generate_hollow_cylinder_mesh"
  "meshcleaner"
  "hex_VS_triangulation_intersection"
  "meshref"
)

cat <<'EOF'
Stage legend:
[0] start
[1] mindist
[2] coarse mesher
[3] mesh filter
[4] monitor function
[5] mesh refinement
[6] fine mesh filter
EOF

if [[ ! -d "${cases_dir}" ]]; then
  echo "Cases directory '${cases_dir}' does not exist" >&2
  exit 1
fi

missing_bins=0
for bin in "${required_binaries[@]}"; do
  bin_path="${script_dir}/${bin}"
  if [[ -x "${bin_path}" ]]; then
    echo "[${bin}][YES]"
  else
    echo "[${bin}][NO]" >&2
    missing_bins=1
  fi
done
if (( missing_bins != 0 )); then
  echo "Required binaries are missing. Ensure all tools are installed next to run_all_cases.sh." >&2
  exit 1
fi

shopt -s nullglob
for case_path in "${cases_dir}"/*; do
  [[ -d "${case_path}" ]] || continue
  case_name=$(basename -- "${case_path}")
  start_time=$(date +%s.%N)
  printf "[%s]:" "${case_name}"
  stage_output=""
  stage_pipe=$(mktemp)
  rm -f "${stage_pipe}"
  mkfifo "${stage_pipe}"
  "${runner}" -s -f "CASES/${case_name}" -n "${num_proc}" "${clean_flag[@]}" > "${stage_pipe}" &
  runner_pid=$!
  while IFS= read -r stage_line; do
    printf "%s" "${stage_line}"
    stage_output+="${stage_line}"
  done < "${stage_pipe}"
  wait "${runner_pid}"
  runner_status=$?
  rm -f "${stage_pipe}"
  if (( runner_status != 0 )); then
    echo
    echo "${stage_output}" >&2
    echo "Case '${case_name}' failed" >&2
    exit 1
  fi
  end_time=$(date +%s.%N)
  elapsed_seconds=$(awk -v start="${start_time}" -v end="${end_time}" 'BEGIN { printf "%.3f", end - start }')
  printf " :: %s [s]\n" "${elapsed_seconds}"
done

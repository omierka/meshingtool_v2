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

num_proc="64"
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
case_paths=()
for case_path in "${cases_dir}"/*; do
  [[ -d "${case_path}" ]] || continue
  case_paths+=("${case_path}")
done

if (( ${#case_paths[@]} == 0 )); then
  echo "No cases found under '${cases_dir}'" >&2
  exit 1
fi

case_names=()
case_name_width=0
for case_path in "${case_paths[@]}"; do
  case_name=$(basename -- "${case_path}")
  case_names+=("${case_name}")
  case_len=${#case_name}
  if (( case_len > case_name_width )); then
    case_name_width=${case_len}
  fi
done

time_placeholder="9999.999 [s]"
time_column_width=${#time_placeholder}
mindist_placeholder="mindist=999.999"
mindist_column_width=${#mindist_placeholder}
nel_placeholder="99,999,999"
nel_column_width=${#nel_placeholder}
iteration_pad_reference="*[2][3][4][5][6]"
iteration_pad_len=${#iteration_pad_reference}

for idx in "${!case_paths[@]}"; do
  case_path="${case_paths[$idx]}"
  case_name="${case_names[$idx]}"
  start_time=$(date +%s.%N)
  printf "[%-${case_name_width}s]:" "${case_name}"
  stage_output=""
  span_value=""
  stage_pipe=$(mktemp)
  rm -f "${stage_pipe}"
  mkfifo "${stage_pipe}"
  "${runner}" -s -f "CASES/${case_name}" -n "${num_proc}" "${clean_flag[@]}" > "${stage_pipe}" &
  runner_pid=$!
  mindist_value=""
  saw_second_iteration=0
  while IFS= read -r stage_line; do
    if [[ "${stage_line}" == mindist=* ]]; then
      mindist_value="${stage_line#mindist=}"
      stage_output+="${stage_line}"$'\n'
      continue
    elif [[ "${stage_line}" == histogram_span=* ]]; then
      span_value="${stage_line#histogram_span=}"
      stage_output+="${stage_line}"$'\n'
      continue
    elif [[ "${stage_line}" == CoarseMeshSize=* ]]; then
      stage_output+="${stage_line}"$'\n'
      continue
    elif [[ "${stage_line}" == "*" ]]; then
      saw_second_iteration=1
    fi
    printf "%s" "${stage_line}"
    stage_output+="${stage_line}"$'\n'
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

  if (( ${#clean_flag[@]} > 0 )); then
    echo
    find "${case_path}" -maxdepth 1 -type f -name '*.png' -delete
    rm -f "${case_path}/size_distribution_histogram.txt"
    continue
  fi

  nel_display=""
  monitor_summary=""
  if (( ${#clean_flag[@]} == 0 )); then
    filtered_tri="${case_path}/Filtered.tri"
    if [[ ! -f "${filtered_tri}" ]]; then
      echo
      echo "Filtered mesh '${filtered_tri}' missing, cannot read NEL/NVT" >&2
      exit 1
    fi
    mesh_line=$(grep 'NEL,NVT' "${filtered_tri}" || true)
    if [[ -z "${mesh_line}" ]]; then
      echo
      echo "Unable to parse NEL/NVT markings in '${filtered_tri}'" >&2
      exit 1
    fi
    # Extract first two integer columns as NEL and NVT counts
    read -r nel nvt _ <<<"${mesh_line}"
    read -r nel_fmt nvt_fmt < <(
      python - "$nel" "$nvt" <<'PY'
import sys
values = [int(v) for v in sys.argv[1:]]
print(" ".join(f"{v:,}" for v in values))
PY
    )
    nel_display=$(printf "NEL=%*s NVT=%*s" "${nel_column_width}" "${nel_fmt}" "${nel_column_width}" "${nvt_fmt}")
    summary_file="${case_path}/monitor_summary.txt"
    if [[ -f "${summary_file}" ]]; then
      monitor_summary=$(<"${summary_file}")
    fi
  fi
  end_time=$(date +%s.%N)
  elapsed_seconds=$(awk -v start="${start_time}" -v end="${end_time}" 'BEGIN { printf "%.3f", end - start }')
  time_str=$(printf "%s [s]" "${elapsed_seconds}")
  printf " :: "
  if (( saw_second_iteration == 0 )); then
    printf "%*s" "${iteration_pad_len}" ""
  fi
  printf "%${time_column_width}s" "${time_str}"
  if (( ${#clean_flag[@]} == 0 )) && [[ -n "${mindist_value}" ]]; then
    mindist_formatted=$(awk -v val="${mindist_value}" 'BEGIN { printf "%.3f", val }')
    if [[ -n "${span_value}" ]]; then
      printf "  mindist=%-${mindist_column_width}s span=%s" "${mindist_formatted}" "${span_value}"
    else
      printf "  %-${mindist_column_width}s" "mindist=${mindist_formatted}"
    fi
  fi
  if [[ -n "${nel_display}" ]]; then
    printf "  :: %s" "${nel_display}"
    if [[ -n "${monitor_summary}" ]]; then
      printf "  %s" "${monitor_summary}"
    fi
  fi
  printf "\n"
done

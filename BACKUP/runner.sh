#!/bin/sh -l

# ---- args ----
FOLDER="PROFEX"
CLEANUP=0
NumProc=4
SILENT=0

while [ $# -gt 0 ]; do
  case "$1" in
    -f|--folder)
      shift
      [ $# -gt 0 ] || { echo "Error: -f|--folder requires a value" >&2; exit 2; }
      FOLDER="$1"
      ;;
    -n|--num-proc)
      shift
      [ $# -gt 0 ] || { echo "Error: -n|--num-proc requires a value" >&2; exit 2; }
      NumProc="$1"
      ;;

    -c|--cleanup)
      CLEANUP=1
      ;;
    -s|--silent)
      SILENT=1
      ;;
    -h|--help)
      echo "Usage: $0 [-f|--folder <folder>] [-n|--num-proc <n>] [-c|--cleanup] [-s|--silent]"
      exit 0
      ;;
    *)
      echo "Error: Unknown argument: $1" >&2
      echo "Usage: $0 [-f|--folder <folder>] [-n|--num-proc <n>] [-c|--cleanup] [-s|--silent]"
      exit 2
      ;;
  esac
  shift
done

# Ensure downstream tools can locate the unified preprocessing defaults
RUNNER_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
if [ -z "${PREPROCESSOR_CONFIG:-}" ]; then
  export PREPROCESSOR_CONFIG="${RUNNER_DIR}/preprocessor.cfg"
fi

# Track location of monitor histogram summary
MONITOR_SUMMARY_FILE="${FOLDER}/monitor_summary.txt"
MONITOR_SUMMARY_VOL_FILE="${FOLDER}/monitor_summary_volumetric.txt"

# ---- optional cleanup ----
if [ "$CLEANUP" -eq 1 ]; then
  echo "Cleanup enabled. Removing generated files under/for: ${FOLDER}"
  rm -rf \
    "${FOLDER}"/*vtu \
    "${FOLDER}"/*tri \
    "${FOLDER}"/*meshDir* \
    "${FOLDER}"/area.txt \
    "${MONITOR_SUMMARY_FILE}" \
    "${MONITOR_SUMMARY_VOL_FILE}"
    exit 0
fi

# Remove stale histogram output so only the final filter exports it
rm -f "${MONITOR_SUMMARY_FILE}" "${MONITOR_SUMMARY_VOL_FILE}"

exec 3>&1
print_stage() {
  if [ "$SILENT" -eq 1 ]; then
    printf '[%s]\n' "$1" >&3
  fi
}

emit_info() {
  if [ "$SILENT" -eq 1 ]; then
    printf '%s\n' "$1" >&3
  else
    printf '%s\n' "$1"
  fi
}

print_iteration_marker() {
  if [ "$SILENT" -eq 1 ]; then
    printf '*\n' >&3
  fi
}

run_mesh_workflow() {
  local coarse_mesh_size="$1"
  local current_span="$2"
  local mindist_value="$3"

  mkdir -p "${FOLDER}/Coarse_meshDir"

  ./generate_hollow_cylinder_mesh \
    -i "${FOLDER}/setup.e3d" \
    -o "${FOLDER}/Coarse_meshDir/Mesh.tri" \
    -s "${coarse_mesh_size}"

  #stage[2]
  print_stage 2

  mpirun -np ${NumProc} ./meshcleaner \
    -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
    -t "${FOLDER}/surface.off" \
    -s 1.0  \
    -o "${FOLDER}"

  mv "${FOLDER}/Coarse_meshDir/Mesh.tri" "${FOLDER}/Coarse_meshDir/Mesh_BU.tri"
  mv "${FOLDER}/Filtered.tri" "${FOLDER}/Coarse_meshDir/Mesh.tri"

  #stage[3]
  print_stage 3

  mpirun -np ${NumProc} ./hex_VS_triangulation_intersection \
   -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
   -t "${FOLDER}/MINGAP.vtu" \
   -m "${mindist_value}" \
   -o "${FOLDER}"

  echo "1" > mesh_names.offs
  echo "${FOLDER}/surface.off" >> mesh_names.offs

  #stage[4]
  print_stage 4

  ./meshref -f "${FOLDER}" \
	    -d "${current_span}"

  #stage[5]
  print_stage 5

  local ispan
  ispan="$(python - "$current_span" <<'PY'
import sys
span = int(float(sys.argv[1]))
print(max(span - 1, 0))
PY
)"

  echo "${FOLDER}/RefinedCleanMesh_lvl${ispan}_refined_clean.vtu"

  rm -f "${MONITOR_SUMMARY_FILE}" "${MONITOR_SUMMARY_VOL_FILE}"
  mpirun -np ${NumProc} ./meshcleaner \
    -h "${FOLDER}/RefinedCleanMesh_lvl${ispan}_refined_clean.vtu" \
    -t "${FOLDER}/surface.off" \
    -s 1.0 \
    -o "${FOLDER}"

  #stage[6]
  print_stage 6
}

read_monitor_zero_percent() {
  local summary_file="$1"
  python - "$summary_file" <<'PY'
import os, re, sys
path = sys.argv[1]
if not os.path.exists(path):
    print("")
    sys.exit(0)
with open(path, "r", encoding="utf-8") as fh:
    data = fh.read()
match = re.search(r'0=\[\s*([0-9.]+)%\]', data)
print(match.group(1) if match else "")
PY
}

maybe_apply_monitor_correction() {
  local count_zero volum_zero need_correction corrected_coarse corrected_span
  count_zero="$(read_monitor_zero_percent "${MONITOR_SUMMARY_FILE}")"
  volum_zero="$(read_monitor_zero_percent "${MONITOR_SUMMARY_VOL_FILE}")"

  if [ -z "${count_zero}" ] || [ -z "${volum_zero}" ]; then
    return
  fi

  need_correction="$(python - "$count_zero" "$volum_zero" <<'PY'
import sys
vals = [float(v) for v in sys.argv[1:]]
eps = 1e-12
print("1" if all(abs(v) <= eps for v in vals) else "0")
PY
)"

  if [ "${need_correction}" != "1" ]; then
    return
  fi

  corrected_coarse="$(python - "$CoarseMeshSize" <<'PY'
import sys
value = float(sys.argv[1]) / 3.0
print(f"{value:.12g}")
PY
)"

  corrected_span="$(python - "$span" <<'PY'
import sys
span = int(float(sys.argv[1]))
print(max(span - 1, 1))
PY
)"

  CoarseMeshSize="${corrected_coarse}"
  span="${corrected_span}"

  emit_info "histogram_span=${span}"

  print_iteration_marker

  run_mesh_workflow "${CoarseMeshSize}" "${span}" "${mindist}"
}

if [ "$SILENT" -eq 1 ]; then
  exec >/dev/null
fi

# ---- environment ----
module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5

#stage [0]
print_stage 0

# ---- compute mindist + CoarseMeshSize ----
read -r mindist span CoarseMeshSize <<EOF
$(./meshhexer-cli --checkpoint-path "${FOLDER}/MINGAP.vtu" min-gap "${FOLDER}/surface.off")
EOF
emit_info "mindist=${mindist}"
emit_info "histogram_span=${span}"
emit_info "CoarseMeshSize=${CoarseMeshSize}"

#stage[1]
print_stage 1

run_mesh_workflow "${CoarseMeshSize}" "${span}" "${mindist}"
maybe_apply_monitor_correction

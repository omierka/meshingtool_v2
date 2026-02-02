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

# ---- optional cleanup ----
if [ "$CLEANUP" -eq 1 ]; then
  echo "Cleanup enabled. Removing generated files under/for: ${FOLDER}"
  rm -rf \
    "${FOLDER}"/*vtu \
    "${FOLDER}"/*tri \
    "${FOLDER}"/*meshDir* \
    "${FOLDER}"/area.txt
    exit 0
fi

exec 3>&1
print_stage() {
  if [ "$SILENT" -eq 1 ]; then
    printf '[%s]\n' "$1" >&3
  fi
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
mindist="$(./meshhexer-cli --checkpoint-path "${FOLDER}/MINGAP.vtu" min-gap "${FOLDER}/surface.off")"
echo "mindist=${mindist}"

#stage[1]
print_stage 1

CoarseMeshSize="$(python -c "print(float('${mindist}')*16.0)")"
echo "CoarseMeshSize=${CoarseMeshSize}"

# ---- workflow ----
mkdir -p "${FOLDER}/Coarse_meshDir"

# Coarse-coarse-mesh-creation
./generate_hollow_cylinder_mesh \
  -i "${FOLDER}/setup.e3d" \
  -o "${FOLDER}/Coarse_meshDir/Mesh.tri" \
  -s "${CoarseMeshSize}"

#stage[2]
print_stage 2
# Mesh-Filtering
mpirun -np ${NumProc} ./meshcleaner \
  -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
  -t "${FOLDER}/surface.off" \
  -s 1.0 \
  -o "${FOLDER}"


mv "${FOLDER}/Coarse_meshDir/Mesh.tri" "${FOLDER}/Coarse_meshDir/Mesh_BU.tri"
mv "${FOLDER}/Filtered.tri" "${FOLDER}/Coarse_meshDir/Mesh.tri"

#stage[3]
print_stage 3
# Monitorfunction-creation
mpirun -np ${NumProc} ./hex_VS_triangulation_intersection \
 -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
 -t "${FOLDER}/MINGAP.vtu" \
 -m "${mindist}" \
 -o "${FOLDER}"

echo "1" > mesh_names.offs
echo "${FOLDER}/surface.off" >> mesh_names.offs

#stage[4]
print_stage 4
# MeshRefinement
./meshref -f "${FOLDER}"

#stage[5]
print_stage 5
# Final-Mesh-Filtering
mpirun -np ${NumProc} ./meshcleaner \
  -h "${FOLDER}/meshDir_BU/Merged_Mesh.tri" \
  -t "${FOLDER}/surface.off" \
  -s 10.0 \
  -o "${FOLDER}"

#stage[6]
print_stage 6

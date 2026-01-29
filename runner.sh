#!/bin/sh -l

# ---- args ----
FOLDER="PROFEX"
CLEANUP=0

while [ $# -gt 0 ]; do
  case "$1" in
    -f|--folder)
      shift
      [ $# -gt 0 ] || { echo "Error: -f|--folder requires a value" >&2; exit 2; }
      FOLDER="$1"
      ;;
    -c|--cleanup)
      CLEANUP=1
      ;;
    -h|--help)
      echo "Usage: $0 [-f|--folder <folder>] [-c|--cleanup]"
      exit 0
      ;;
    *)
      echo "Error: Unknown argument: $1" >&2
      echo "Usage: $0 [-f|--folder <folder>] [-c|--cleanup]"
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

# ---- environment ----
module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5

# ---- compute mindist + CoarseMeshSize ----
mindist="$(meshhexer-cli --checkpoint-path "${FOLDER}/MINGAP.vtu" min-gap "${FOLDER}/surface.off")"
echo "mindist=${mindist}"

CoarseMeshSize="$(python -c "print(float('${mindist}')*16.0)")"
echo "CoarseMeshSize=${CoarseMeshSize}"

# ---- workflow ----
mkdir -p "${FOLDER}/Coarse_meshDir"

# Coarse-coarse-mesh-creation
generate_hollow_cylinder_mesh.py \
  -i "${FOLDER}/MESHCONFIG.dat" \
  -o "${FOLDER}/Coarse_meshDir/Mesh.tri" \
  -s "${CoarseMeshSize}"

# Mesh-Filtering
meshcleaner \
  -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
  -t "${FOLDER}/surface.off" \
  -s 1.0 \
  -o "${FOLDER}"


mv "${FOLDER}/Coarse_meshDir/Mesh.tri" "${FOLDER}/Coarse_meshDir/Mesh_BU.tri"
mv "${FOLDER}/Filtered.tri" "${FOLDER}/Coarse_meshDir/Mesh.tri"

# Monitorfunction-creation
mpirun -np 64 hex_VS_triangulation_intersection \
 -h "${FOLDER}/Coarse_meshDir/Mesh.tri" \
 -t "${FOLDER}/MINGAP.vtu" \
 -m "${mindist}" \
 -o "${FOLDER}"

echo "1" > mesh_names.offs
echo "${FOLDER}/surface.off" >> mesh_names.offs

# MeshRefinement
meshref -f "${FOLDER}"

# Final-Mesh-Filtering
meshcleaner \
  -h "${FOLDER}/NEW_meshDir/ReducedMesh.tri" \
  -t "${FOLDER}/surface.off" \
  -s 10.0 \
  -o "${FOLDER}"

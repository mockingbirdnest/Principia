# This script searches for all copies of KSP on your MacOS computer and adds
# references to them in ../KSP Assemblies. This will allow you to build the C#
# project.

# The script utilizes MacOS's metadata utilities and its application structure.
# It does not work on other OSes.
AGENT_OS=$(uname -s)
if [[ "${AGENT_OS}" != "Darwin" ]]; then
  echo "This script only supports MacOS."
  exit 1
fi

# Create the KSP Assemblies dir if not present.
BASE_DIR=$(dirname $0)
KSP_ASSEMBLIES_DIR="${BASE_DIR}/../KSP Assemblies"
mkdir -p "${KSP_ASSEMBLIES_DIR}"

# Get all installed copies of KSP.
KSP_BUNDLE_IDENTIFIER="unity.Squad.Kerbal Space Program"
IFS=$'\n' KSP_APPS=($(mdfind "kMDItemCFBundleIdentifier == '${KSP_BUNDLE_IDENTIFIER}'"))

for KSP_APP in "${KSP_APPS[@]}"; do
  echo "Found ${KSP_APP}"
  
  # Get the KSP version.
  KSP_VERSION=$(mdls -name=kMDItemVersion -raw "${KSP_APP}") 
  echo "Version: ${KSP_VERSION}"

  # Link the KSP app's dlls to the KSP Assemblies dir.
  SOURCE_DIR="${KSP_APP}/Contents/Resources/Data/Managed"
  TARGET_DIR="${KSP_ASSEMBLIES_DIR}/${KSP_VERSION}"
  ln -s "${SOURCE_DIR}" "${TARGET_DIR}"
done

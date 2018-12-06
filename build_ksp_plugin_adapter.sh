#!/bin/bash
# Please specify the KSP_MANAGED_FOLDER variable if the location is different from the default in the script
# Typically this is something like: "${LOCATION_OF_GAMES}/Kerbal Space Program/KSP_Data/Managed/"
# This folder contains the KSP and unity C# assembly to build against

if [ -z ${KSP_MANAGED_FOLDER} ]; then
	KSP_MANAGED_FOLDER="${HOME}/.steam/steam/steamapps/common/Kerbal Space Program/KSP_Data/Managed/"
fi

AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" msbuild ksp_plugin_adapter/ksp_plugin_adapter.csproj /t:Build /p:Configuration="Debug"
AssemblySearchPaths="/usr/lib/mono/3.5-api/;/usr/lib/mono/2.0-api/;${KSP_MANAGED_FOLDER}" msbuild ksp_plugin_adapter/ksp_plugin_adapter.csproj /t:Build /p:Configuration="Release"

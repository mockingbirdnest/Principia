# It seems that protoc really wants its dependencies to be in /usr/local/lib.
# In some setups, e.g., Azure pipelines, this does not work, so we need to help
# it find its dynamic libraries.
if [ "${AGENT_OS?}" == "Darwin" ]; then
  install_name_tool \
      -change \
          /usr/local/lib/libprotobuf.17.dylib \
          ./deps/protobuf/src/.libs/libprotobuf.17.dylib \
      deps/protobuf/src/.libs/libprotoc.17.dylib
  install_name_tool \
      -change \
          /usr/local/lib/libprotobuf.17.dylib \
          ./deps/protobuf/src/.libs/libprotobuf.17.dylib \
      deps/protobuf/src/.libs/protoc
  install_name_tool \
      -change \
          /usr/local/lib/libprotoc.17.dylib \
          ./deps/protobuf/src/.libs/libprotoc.17.dylib \
      deps/protobuf/src/.libs/protoc
elif [ "${AGENT_OS?}" == "Linux" ]; then
  export LD_LIBRARY_PATH="./deps/protobuf/src/.libs:$LD_LIBRARY_PATH"
fi

make test
make release

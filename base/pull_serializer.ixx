module;

#include "absl/base/thread_annotations.h"
#include "absl/synchronization/mutex.h"
#include "gipfeli/compression.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/message.h"
// For body:
#include "absl/log/check.h"

export module principia.base.pull_serializer;

#include "base/pull_serializer.hpp"

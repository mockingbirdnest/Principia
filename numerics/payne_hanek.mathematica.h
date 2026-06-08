#pragma once

#include <array>
#include <cstdint>

namespace principia {
namespace numerics {
namespace _payne_hanek {
namespace internal {

constexpr std::int64_t PayneHanekBitsPerChunk = 26;

// Chunks of 26 bits of 4/π, up to exponent -1074.
constexpr std::array<double, 42> PayneHanekChunks{
    0x28b'e60d.0p-25,
    0x2e4'e441.0p-51,
    0x14a'7f09.0p-77,
    0x357'd1f5.0p-103,
    0x0d3'7703.0p-129,
    0x1b6'2959.0p-155,
    0x24f'10e4.0p-181,
    0x041'fe51.0p-207,
    0x18e'af7a.0p-233,
    0x3bc'561b.0p-259,
    0x1c9'1b8e.0p-285,
    0x242'4dd2.0p-311,
    0x380'1924.0p-337,
    0x2ee'a09d.0p-363,
    0x064'873f.0p-389,
    0x21d'eb1c.0p-415,
    0x2c4'a69c.0p-441,
    0x3ee'8823.0p-467,
    0x17d'4bae.0p-493,
    0x344'84e9.0p-519,
    0x271'c09a.0p-545,
    0x345'f7e4.0p-571,
    0x04e'6475.0p-597,
    0x239'8353.0p-623,
    0x0e7'd272.0p-649,
    0x045'f8bb.0p-675,
    0x37e'4a0e.0p-701,
    0x31f'f897.0p-727,
    0x3ff'7816.0p-753,
    0x180'fef2.0p-779,
    0x3c4'62d6.0p-805,
    0x20a'6d1f.0p-831,
    0x1b4'd9fb.0p-857,
    0x0f2'7cb0.0p-883,
    0x26d'd3d1.0p-909,
    0x23f'669e.0p-935,
    0x17f'a8b5.0p-961,
    0x352'7bac.0p-987,
    0x1fa'f97c.0p-1013,
    0x17b'3d07.0p-1039,
    0x0e7'de29.0p-1065,
    0x129'2ea6.0p-1091};

}  // namespace internal

using internal::PayneHanekBitsPerChunk;
using internal::PayneHanekChunks;

}  // namespace _payne_hanek
}  // namespace numerics
}  // namespace principia

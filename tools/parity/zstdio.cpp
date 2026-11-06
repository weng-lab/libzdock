// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Minimal zstd stream filter for the parity snapshot tooling: compresses
// or decompresses stdin to stdout (zstd format), linked against the
// system libzstd so snapshot.py needs no zstd CLI or Python module.
//
//   zstdio c  < input > output.zst
//   zstdio d  < input.zst > output

#include <zstd.h>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace {

constexpr std::size_t kBufSize = 1 << 20;

void fail(const char *what) {
  std::fprintf(stderr, "zstdio: %s\n", what);
  std::exit(1);
}

void failZstd(const char *what, std::size_t code) {
  std::fprintf(stderr, "zstdio: %s: %s\n", what, ZSTD_getErrorName(code));
  std::exit(1);
}

void run(bool compress) {
  static std::string ibuf(kBufSize, '\0');
  static std::string obuf(kBufSize * 4, '\0');
  if (compress) {
    ZSTD_CCtx *ctx = ZSTD_createCCtx();
    if (!ctx)
      fail("cannot create compression context");
    ZSTD_CCtx_setParameter(ctx, ZSTD_c_compressionLevel, 3);
    while (true) {
      const std::size_t got = std::fread(ibuf.data(), 1, ibuf.size(), stdin);
      if (std::ferror(stdin))
        fail("read error");
      const bool last = std::feof(stdin) != 0;
      ZSTD_inBuffer in{ibuf.data(), got, 0};
      while (true) {
        ZSTD_outBuffer out{obuf.data(), obuf.size(), 0};
        const std::size_t remaining = ZSTD_compressStream2(
            ctx, &out, &in, last ? ZSTD_e_end : ZSTD_e_continue);
        if (ZSTD_isError(remaining))
          failZstd("compress", remaining);
        if (out.pos > 0 &&
            std::fwrite(obuf.data(), 1, out.pos, stdout) != out.pos)
          fail("write error");
        if (last && remaining == 0)
          break; // everything in and flushed
        if (!last && in.pos >= in.size)
          break; // this chunk consumed
      }
      if (last)
        break;
    }
    ZSTD_freeCCtx(ctx);
  } else {
    ZSTD_DCtx *ctx = ZSTD_createDCtx();
    if (!ctx)
      fail("cannot create decompression context");
    while (true) {
      const std::size_t got = std::fread(ibuf.data(), 1, ibuf.size(), stdin);
      if (got == 0)
        break;
      if (std::ferror(stdin))
        fail("read error");
      ZSTD_inBuffer in{ibuf.data(), got, 0};
      while (in.pos < in.size) {
        ZSTD_outBuffer out{obuf.data(), obuf.size(), 0};
        const std::size_t ret = ZSTD_decompressStream(ctx, &out, &in);
        if (ZSTD_isError(ret))
          failZstd("decompress", ret);
        if (out.pos > 0 &&
            std::fwrite(obuf.data(), 1, out.pos, stdout) != out.pos)
          fail("write error");
      }
    }
    ZSTD_freeDCtx(ctx);
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2 || (std::string(argv[1]) != "c" &&
                    std::string(argv[1]) != "d")) {
    std::fprintf(stderr, "usage: zstdio c|d  (stdin -> stdout)\n");
    return 2;
  }
  run(argv[1][0] == 'c');
  return 0;
}

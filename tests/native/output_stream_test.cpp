// Adapted from RAD dev 6167540 for the embedded standard-queue/gzip backend.
#include "rad/rad_headers.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <system_error>
#include <sys/stat.h>
#include <unistd.h>

using namespace rrad_rad_core;

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

std::string contents(FILE* file) {
    require(std::fseek(file, 0, SEEK_SET) == 0, "caller-owned output stream was closed");
    std::string result;
    char buffer[256];
    for (size_t n; (n = std::fread(buffer, 1, sizeof(buffer), file)) != 0;)
        result.append(buffer, n);
    require(!std::ferror(file), "cannot read output stream fixture");
    return result;
}

void small_and_empty_streams() {
    FILE* file = std::tmpfile();
    require(file != nullptr, "cannot create output stream fixture");
    const std::string record = "@tiny\tCB:Z:ACGT\nACGT\n+\nIJKL\n";
    {
        sigstring_writing sink(file, sigstring_writing::format::FASTQA);
        {
            parallel_writer writer;
            writer.write_raw_string(sink, std::string(record));
            writer.stop();
        }
        sink.finish();
        sink.finish();
        require(contents(file) == record, "small FASTQ output or final tail was lost");
    }
    require(contents(file) == record, "writer destructor closed or changed a non-owning stream");
    require(std::fclose(file) == 0, "caller could not close output stream");

    file = std::tmpfile();
    require(file != nullptr, "cannot create empty output stream fixture");
    {
        sigstring_writing sink(file, sigstring_writing::format::FASTQA);
        sink(std::string{});
        sink.finish();
    }
    require(contents(file).empty(), "empty stream emitted non-record bytes");
    require(std::fclose(file) == 0, "cannot close empty output stream");
}

void asynchronous_write_failure() {
    FILE* file = std::fopen("/dev/full", "wb");
    if (!file) {
#ifdef __linux__
        throw std::runtime_error("Linux output-failure fixture /dev/full is unavailable");
#else
        std::cout << "Platform fixture: /dev/full unavailable; Linux CI covers ENOSPC output failure\n";
        return;
#endif
    }
    bool saw_failure = false;
    {
        sigstring_writing sink(file, sigstring_writing::format::FASTQA);
        parallel_writer writer;
        try {
            // More than the queue capacity also exercises failed-producer
            // cancellation while a slow/failed sink cannot accept more jobs.
            for (int i = 0; i < 256; ++i)
                writer.write_raw_string(sink, std::string(65536, 'A'));
            writer.stop();
        } catch (const std::system_error& error) {
            require(error.code().value() == ENOSPC || error.code().value() == EIO,
                    "writer lost the output failure's error code");
            saw_failure = true;
        }
        try { writer.stop(); }
        catch (const std::system_error&) { saw_failure = true; }
        try { sink.finish(); }
        catch (const std::system_error&) { saw_failure = true; }
    }
    std::fclose(file);
    require(saw_failure, "asynchronous output failure was silently ignored");
}

void chunk_failure_reaches_caller() {
    char path[] = "/tmp/rad-output-chunks-XXXXXX.fastq";
    const int fd = mkstemps(path, 6);
    require(fd >= 0, "cannot create chunk fixture");
    FILE* file = fdopen(fd, "w");
    require(file != nullptr, "cannot open chunk fixture");
    for (int i = 0; i < 8; ++i)
        require(std::fprintf(file, "@r%d\nACGT\n+\nIIII\n", i) > 0, "cannot write chunk fixture");
    require(std::fclose(file) == 0, "cannot close chunk fixture");
    bool caught = false;
    auto fail = [](std::vector<read_streaming::sequence>&, const std::string&) {
        throw std::runtime_error("expected chunk output failure");
    };
    try {
        chunk_streaming<read_streaming::sequence, decltype(fail)> reader(1, 1);
        reader.process_chunks(path, fail, 2);
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "expected chunk output failure";
    }
    unlink(path);
    require(caught, "chunk-worker exception did not reach the calling thread");
}

void compressor_finish_failure() {
    char directory[] = "/tmp/rad-output-compressor-XXXXXX";
    require(mkdtemp(directory) != nullptr, "cannot create compressor fixture directory");
    const std::string script = std::string(directory) + "/pigz-failure";
    const std::string output = std::string(directory) + "/output.gz";
    {
        std::ofstream fixture(script);
        fixture << "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then exit 0; fi\n"
                   "cat >/dev/null\nexit 7\n";
        require(bool(fixture), "cannot write compressor fixture");
    }
    require(chmod(script.c_str(), 0700) == 0, "cannot enable compressor fixture");
    // No earlier case asks for compression, so the resolver's first cached
    // choice is this executable. The failed exit occurs only after EOF.
    unsetenv("RAD_NO_PIGZ");
    setenv("RAD_PIGZ", script.c_str(), 1);
    bool caught = false;
    try {
        sigstring_writing sink(output, sigstring_writing::format::FASTQA, true, false, 1);
        sink(std::string("@tail\nACGT\n+\nIIII\n"));
        sink.finish();
    } catch (const std::exception&) {
        caught = true;
    }
    unlink(script.c_str());
    unlink(output.c_str());
    rmdir(directory);
    require(caught, "compressor failure on finalization was silently ignored");
}
}

int main() {
    try {
        small_and_empty_streams();
        asynchronous_write_failure();
        chunk_failure_reaches_caller();
        compressor_finish_failure();
        std::cout << "Output streams preserve tails and ownership; asynchronous and chunk errors propagate\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

#pragma once
#include "rad_headers.h"

// needed a way for IO to work on .gz and non .gz files, and needed a PIGZ/non-PIGZ solution for files.
// finangled a way to deal with this by just opening two kseq instantiations because you can't trade them off once opened.
// so just open two, kseq_fd and kseq_gz, and use the appropriate one based on whether pigz is being used or not/filetype.
// this is all wrapped up in file_streaming and file_pointer classes below. Learned a lot about reader threads and async IO.
// basically, the smart way to do things is to have pigz do the decompression in parallel with multiple threads (`pigz_reading`) 
// and then have a reader class (`file_streaming`) that manages the file pointers and kseq instances 
// for reading stuff through kseq.

/**
 * @brief kseq intialization for regular file
 */
namespace kseq_fd {
    #ifndef KSEQ_INIT_INT_READ
    #define KSEQ_INIT_INT_READ
    KSEQ_INIT(int, read)
    #endif
}

/**
 * @brief kseq initialization for gzFile
 */
namespace kseq_gz {
    #ifndef KSEQ_INIT_GZFILE_GZREAD  
    #define KSEQ_INIT_GZFILE_GZREAD
    KSEQ_INIT(gzFile, gzread)
    #endif
}

/**
 * @brief Parallel decompression reading using pigz
 * 
 * Manages pigz subprocess for multi-threaded decompression of gzipped files
 */
class pigz_reading {
public:
    static std::string shell_quote_posix(const std::string& s) {
        // POSIX /bin/sh-safe single-quote escaping:
        // abc'd -> 'abc'"'"'d'
        std::string out;
        out.reserve(s.size() + 2);
        out.push_back('\'');
        for (char c : s) {
            if (c == '\'') {
                out.append("'\"'\"'");
            } else {
                out.push_back(c);
            }
        }
        out.push_back('\'');
        return out;
    }

    static std::string sanitize_user_path(std::string p) {
        // Trim ASCII whitespace (incl. CR/LF) commonly introduced by CSV/TSV parsing.
        const auto start = p.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        const auto end = p.find_last_not_of(" \t\r\n");
        p = p.substr(start, end - start + 1);

        // Strip UTF-8 BOM if present.
        if (p.size() >= 3) {
            const unsigned char b0 = static_cast<unsigned char>(p[0]);
            const unsigned char b1 = static_cast<unsigned char>(p[1]);
            const unsigned char b2 = static_cast<unsigned char>(p[2]);
            if (b0 == 0xEF && b1 == 0xBB && b2 == 0xBF) {
                p.erase(0, 3);
            }
        }

        // Strip surrounding quotes if the whole field is quoted.
        if (p.size() >= 2) {
            const char f = p.front();
            const char l = p.back();
            if ((f == '"' && l == '"') || (f == '\'' && l == '\'')) {
                p = p.substr(1, p.size() - 2);
            }
        }

        // Expand leading "~/" to $HOME (shell would normally do this).
        if (!p.empty() && p[0] == '~' && (p.size() == 1 || p[1] == '/')) {
            if (const char* home = std::getenv("HOME"); home && *home) {
                p = std::string(home) + p.substr(1);
            }
        }
        return p;
    }

    /**
     * @brief Open a pigz decompression pipe for reading
     * @param input_path Path to compressed file (.gz, .fq.gz, .fastq.gz, etc.)
     * @param threads Number of decompression threads (default: 4)
     * @return FILE* pointer to pipe, or nullptr on failure
     */
    inline FILE* open_pigz_read_pipe(const std::string& input_path, int threads = 4) {
        // Resolve pigz (PATH-first; honors RAD_PIGZ / RAD_NO_PIGZ). "" => skip.
        const std::string& pigz = resolve_pigz();
        if (pigz.empty()) {
            return nullptr;
        }

        const std::string path = sanitize_user_path(input_path);
        if (path.empty()) {
            return nullptr;
        }

        // Avoid invoking pigz on a path that doesn't exist; this also makes it
        // much easier to diagnose hidden characters (e.g. CR, BOM) in the path.
        std::error_code ec;
        if (!std::filesystem::is_regular_file(path, ec)) {
            return nullptr;
        }

        // Use pigz -dc for decompression to stdout
        // -d = decompress, -c = write to stdout, -p = threads
        std::string cmd = shell_quote_posix(pigz) + " -dc -p " + std::to_string(threads) +
                          " " + shell_quote_posix(path);
        
        FILE* fp = popen(cmd.c_str(), "r");
        if (!fp) {
            return nullptr;
        }
        
        //removed because if you set a static thread_local buffer, it messes with multiple files opened in different threads.
        // will cause issues streaming from multiple files in parallel
        // had originally set it up trying to debug, TSPMTFO

        // Large read buffer (16 MB) to reduce read() syscalls
        
        //static thread_local std::vector<char> read_buf(16u << 20);
       // setvbuf(fp, read_buf.data(), _IOFBF, read_buf.size());
        
        return fp;
    }

    /**
     * @brief Close pigz decompression pipe
     * @param fp FILE pointer from open_pigz_read_pipe
     * @return Exit status of pigz process
     */
    inline int close_pigz_read_pipe(FILE* fp) {
        if (!fp) return -1;
        return pclose(fp);
    }

    /**
     * @brief Resolve the pigz executable to use, preferring one on PATH (so a
     * conda or system pigz "just works"), then the compile-time PIGZ_PATH build
     * hint. Honors RAD_PIGZ (explicit override) and RAD_NO_PIGZ (disable).
     * Returns "" when pigz is unavailable or disabled. Resolved once + cached.
     */
    static const std::string& resolve_pigz() {
        static const std::string cached = []() -> std::string {
            if (const char* no = std::getenv("RAD_NO_PIGZ"); no && *no) return "";
            auto available = [](const std::string& exe) -> bool {
                if (exe.empty()) return false;
                std::string c = "command -v " + shell_quote_posix(exe) + " > /dev/null 2>&1";
                return std::system(c.c_str()) == 0;
            };
            // explicit override wins (if it resolves)
            if (const char* env = std::getenv("RAD_PIGZ"); env && *env) {
                return available(env) ? std::string(env) : std::string();
            }
            // pigz on PATH -- the common case (conda run-dep, system install)
            if (available("pigz")) return "pigz";
            // fall back to the compile-time build hint if it still resolves
            #ifdef PIGZ_PATH
            { std::string p = PIGZ_PATH; if (p != "pigz" && available(p)) return p; }
            #endif
            return "";
        }();
        return cached;
    }

    /**
     * @brief Check if pigz is available (on PATH, via RAD_PIGZ, or PIGZ_PATH)
     * @return true if a usable pigz was found and is not disabled
     */
    static bool is_pigz_available() {
        return !resolve_pigz().empty();
    }
};

/**
 * @brief File streaming with automatic parallel decompression
 * 
 * Automatically uses pigz for .gz files when available, falls back to gzip.
 * Handles both single files and directories of FASTQ/FASTA files.
 */
class file_streaming {
public:
/**
 * @brief File pointer structure for kseq reading
 * @param fp_pigz Pigz pipe handle or regular file handle
 * @param fp_gzip Fallback gzip handle
 * @param seq_fd kseq for FILE* (pigz path or uncompressed)
 * @param seq_gz kseq for gzFile (fallback path)
 * @param fd File descriptor for kseq
 * @param path File path
 * @param using_pigz Flag indicating if pigz is used
 * @param pigz_reader pigz reading manager
 */
    struct file_pointer {
        FILE* fp_pigz;      // Pigz pipe handle or regular file handle
        gzFile fp_gzip;     // Fallback gzip handle
        kseq_fd::kseq_t* seq_fd;     // kseq for FILE* (pigz path or uncompressed)
        kseq_gz::kseq_t* seq_gz;     // kseq for gzFile (fallback path)
        int fd;             // File descriptor for kseq
        std::string path;
        bool using_pigz;
        bool owns_pigz_pipe;
        pigz_reading pigz_reader;
        
        inline kseq_fd::kseq_t* get_seq_fd() { 
            return seq_fd; 
        }

        inline kseq_gz::kseq_t* get_seq_gz() { 
            return seq_gz; 
        }

        // one-step read that picks the right backend
        inline int kseq_read_any() {
            if (using_pigz) {
                return seq_fd ? kseq_fd::kseq_read(seq_fd) : -1;  // -1 = EOF, -2 = bad FASTQ
            } else {
                return seq_gz ? kseq_gz::kseq_read(seq_gz) : -1;
            }
        }

        // unified field accessors (avoid duplicating branches)
        inline const char* name_s() const {
            return using_pigz ? seq_fd->name.s : seq_gz->name.s;
        }
        inline const char* comment_s() const {
            return using_pigz ? (seq_fd->comment.l ? seq_fd->comment.s : nullptr)
                            : (seq_gz->comment.l ? seq_gz->comment.s : nullptr);
        }
        inline bool is_fastq_record() const {
            return using_pigz ? seq_fd->is_fastq : seq_gz->is_fastq;
        }
        inline const char* plus_s() const {
            return using_pigz ? (seq_fd->plus.l ? seq_fd->plus.s : nullptr)
                            : (seq_gz->plus.l ? seq_gz->plus.s : nullptr);
        }
        inline const char* seq_s() const {
            return using_pigz ? seq_fd->seq.s : seq_gz->seq.s;
        }
        inline int qual_len() const {
            return using_pigz ? seq_fd->qual.l : seq_gz->qual.l;
        }
        inline const char* qual_s() const {
            return using_pigz ? (seq_fd->qual.l ? seq_fd->qual.s : nullptr)
                            : (seq_gz->qual.l ? seq_gz->qual.s : nullptr);
        }

        /**
         * @brief Constructor to open file pointer
         * @param file_path Path to input file
         * @param pigz_threads Number of pigz threads for decompression
         * @note Automatically detects .gz files and uses pigz/gzip as needed
         */
        file_pointer(const std::string& file_path, int pigz_threads = 4,
                     bool allow_pigz = true)
            : fp_pigz(nullptr), fp_gzip(nullptr), seq_fd(nullptr), seq_gz(nullptr),
              fd(-1), path(file_path), using_pigz(false),
              owns_pigz_pipe(false) {
            path = pigz_reading::sanitize_user_path(path);
            if (path.empty()) {
                return;
            }
            
            // Check if file is compressed
            bool is_compressed = has_gz_extension(path);
            
            if (is_compressed) {
                // Try pigz first for parallel decompression
                if (allow_pigz && pigz_reading::is_pigz_available()) {
                    fp_pigz = pigz_reader.open_pigz_read_pipe(path, pigz_threads);
                    if (fp_pigz) {
                        fd = fileno(fp_pigz);
                        seq_fd = kseq_fd::kseq_init(fd);
                        using_pigz = true;
                        owns_pigz_pipe = true;
                        return;
                    }
                }

                // Fallback to single-threaded gzip
                fp_gzip = gzopen(path.c_str(), "r");
                if (fp_gzip) {
                    gzbuffer(fp_gzip, 1 << 22);  // 4 MiB buffer
                    seq_gz = kseq_gz::kseq_init(fp_gzip);
                    using_pigz = false;
                }
            } else {
                // Uncompressed file - use regular file descriptor
                fp_pigz = fopen(path.c_str(), "r");
                if (fp_pigz) {
                    fd = fileno(fp_pigz);
                    seq_fd =  kseq_fd::kseq_init(fd);
                    setvbuf(fp_pigz, nullptr, _IOFBF, 1 << 22);  // 4 MiB buffer
                    using_pigz = true;  // Using FILE*, not gzFile
                    owns_pigz_pipe = false;
                }
            }
        }

        /**
         * @brief Destructor to clean up resources
         */
        ~file_pointer() {
            if (seq_fd) kseq_destroy(seq_fd);
            if (seq_gz) kseq_destroy(seq_gz);
            if (using_pigz && fp_pigz) {
                if (owns_pigz_pipe) {
                    pigz_reader.close_pigz_read_pipe(fp_pigz);
                } else {
                    std::fclose(fp_pigz);
                }
            } else if (fp_gzip) {
                gzclose(fp_gzip);
            }
        }

        bool is_valid() const {
            return (using_pigz && fp_pigz && seq_fd) || (fp_gzip && seq_gz);
        }

        bool gzip_read_failed(int& code, std::string& message) const {
            code = Z_OK;
            message.clear();
            if (using_pigz || !fp_gzip) return false;
            const char* detail = gzerror(fp_gzip, &code);
            if (code == Z_OK || code == Z_STREAM_END) return false;
            if (detail) message = detail;
            return true;
        }

        // kseq reports a clean-looking EOF when the decompressor feeding its
        // pipe exits early (for example, after detecting a bad gzip trailer).
        // Reap pigz at the EOF boundary and turn any non-zero pclose status
        // into a read failure before the caller advances to another file or
        // reports success.  This is deliberately not done by the destructor:
        // callers that stop at max_reads close a still-producing pipe and may
        // legitimately make pigz exit via SIGPIPE.
        bool pigz_read_failed(int& code, std::string& message) {
            code = 0;
            message.clear();
            if (!using_pigz || !owns_pigz_pipe || !fp_pigz) return false;

            if (seq_fd) {
                kseq_destroy(seq_fd);
                seq_fd = nullptr;
            }
            errno = 0;
            FILE* pipe = fp_pigz;
            fp_pigz = nullptr;
            owns_pigz_pipe = false;
            const int status = pigz_reader.close_pigz_read_pipe(pipe);
            const int saved_errno = errno;
            if (status == 0) return false;

            code = status;
            if (status == -1) {
                message = "could not collect pigz exit status";
                if (saved_errno != 0) {
                    message += ": ";
                    message += std::strerror(saved_errno);
                }
            } else {
                message = "pigz decompression exited unsuccessfully (pclose status " +
                          std::to_string(status) + ")";
            }
            return true;
        }

    private:
        static bool has_gz_extension(const std::string& path) {
            if (path.size() < 3) return false;
            const size_t offset = path.size() - 3;
            return path[offset] == '.' &&
                   std::tolower(static_cast<unsigned char>(path[offset + 1])) ==
                       'g' &&
                   std::tolower(static_cast<unsigned char>(path[offset + 2])) ==
                       'z';
        }
    };

    std::unique_ptr<file_pointer> current_file;
    std::filesystem::directory_iterator dir_iter;
    bool dir_status;
    std::string path;
    int pigz_threads;
    bool allow_pigz;

    /**
     * @brief Constructor to initialize file streaming
     * @param input_path Path to input file or directory
     * @param threads Number of pigz threads for decompression
     * @throws std::runtime_error if file cannot be opened or path is invalid
     */
    explicit file_streaming(const std::string& input_path, int threads = 4,
                            bool use_pigz = true)
        : pigz_threads(threads), allow_pigz(use_pigz) {
        path = pigz_reading::sanitize_user_path(input_path);
        dir_status = is_directory(path);
        if (dir_status) {
            dir_iter = std::filesystem::directory_iterator(path);
            open_next_file();
        } else if (is_fastqa(path)) {
            current_file = std::make_unique<file_pointer>(
                path, pigz_threads, allow_pigz);
            if (!current_file->is_valid()) {
                throw std::runtime_error("Failed to open sequence file: " + path);
            }
        } else if (is_file(path)) {
            throw std::invalid_argument(
                "unsupported sequence input: " + path +
                " (expected FASTQ/FASTA, optionally gzip-compressed)");
        } else {
            throw std::invalid_argument("Invalid input path: " + path);
        }
    }

    bool open_next_file() {
        if (!dir_status) return false;
        while (dir_iter != std::filesystem::directory_iterator()) {
            const auto& entry = *dir_iter++;
            if (!is_file(get_file_path(entry))) continue;
            
            std::string file_path = entry.path().string();
            if (is_fastqa(file_path)) {
                current_file = std::make_unique<file_pointer>(
                    file_path, pigz_threads, allow_pigz);
                if (current_file->is_valid()) {
                    return true;
                }
            }
        }
        return false;
    }

    static bool has_suffix(const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() &&
               str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    static bool is_fastqa(const std::string& path) {
        static const std::vector<std::string> extensions = {
            ".fastq", ".fq", ".fasta", ".fa",
            ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz"
        };

        std::string normalized = path;
        std::transform(normalized.begin(), normalized.end(),
                       normalized.begin(), [](unsigned char value) {
                           return static_cast<char>(std::tolower(value));
                       });
        for (const auto& ext : extensions) {
            if (has_suffix(normalized, ext)) return true;
        }
        return false;
    }

    static bool is_directory(const std::string& path) {
        return std::filesystem::is_directory(path);
    }

    static bool is_file(const std::string& path) {
        return std::filesystem::is_regular_file(path);
    }

    static std::string get_file_path(const std::filesystem::directory_entry& entry) {
        return entry.path().string();
    }

    std::string current_file_path() const {
        return current_file ? current_file->path : "";
    }
};

/**
 * @brief Read streaming for FASTQ/FASTA sequences
 */
class read_streaming {
public:
    /**
     * @brief structure of an incoming read
     * @param id original read id `string`
     * @param comment optional comment `string`
     * @param seq  `string` sequence
     * @param plus optional FASTQ plus-line payload, excluding the leading `+`
     * @param qual  `string` qual string for fastq
     * @param is_fastq `bool` flag indicating if read is from FASTQ
     */
    struct sequence {
        std::string id;
        std::string comment;
        std::string seq;
        std::string plus;
        std::string qual;
        bool is_fastq;
    };

    file_streaming& files;
    int read_error_code_ = 0;
    std::string read_error_path_;
    std::string read_error_detail_;
    uint64_t empty_sequences_skipped_ = 0;
    bool preserve_empty_sequences_ = false;
    
    explicit read_streaming(file_streaming& f,
                            bool preserve_empty_sequences = false)
        : files(f), preserve_empty_sequences_(preserve_empty_sequences) {}

    bool read_failed() const {
        return read_error_code_ < -1;
    }

    int read_error_code() const {
        return read_error_code_;
    }

    const std::string& read_error_path() const {
        return read_error_path_;
    }

    const std::string& read_error_detail() const {
        return read_error_detail_;
    }

    uint64_t empty_sequences_skipped() const {
        return empty_sequences_skipped_;
    }

    void note_empty_sequence(file_streaming::file_pointer& fp) {
        ++empty_sequences_skipped_;
        if (empty_sequences_skipped_ <= 10) {
            const char* name = fp.name_s();
            RAD_CERR << "Warning: skipping zero-length sequence record"
                      << (name ? " '" + std::string(name) + "'" : "")
                      << " in " << files.current_file_path() << "\n";
        } else if (empty_sequences_skipped_ == 11) {
            RAD_CERR
                << "Warning: additional zero-length sequence warnings "
                   "suppressed\n";
        }
    }

    /**
     * @brief Get a view of the next sequence without copying read metadata.
     *
     * The view remains valid only until the next read from this object.  This
     * is intended for serial consumers that finish processing a sequence
     * before advancing kseq.
     */
    std::optional<std::string_view> next_sequence_view() {
        while (files.current_file || files.open_next_file()) {
            auto& fp = files.current_file;
            int l = fp->kseq_read_any();
            if (l == 0 && !preserve_empty_sequences_) {
                note_empty_sequence(*fp);
                continue;
            }
            if (l >= 0) {
                return std::string_view(
                    fp->seq_s(), static_cast<size_t>(l));
            }
            int gzip_code = Z_OK;
            std::string gzip_detail;
            if (l < 0 && fp->gzip_read_failed(gzip_code, gzip_detail)) {
                read_error_code_ = l < -1 ? l : -3;
                read_error_path_ = files.current_file_path();
                read_error_detail_ = "gzip/zlib code " +
                                     std::to_string(gzip_code);
                if (!gzip_detail.empty()) {
                    read_error_detail_ += ": " + gzip_detail;
                }
                return std::nullopt;
            }
            int pigz_status = 0;
            std::string pigz_detail;
            if (l < 0 && fp->pigz_read_failed(pigz_status, pigz_detail)) {
                read_error_code_ = -3;
                read_error_path_ = files.current_file_path();
                read_error_detail_ = pigz_detail;
                return std::nullopt;
            }
            if (l < -1) {
                read_error_code_ = l;
                read_error_path_ = files.current_file_path();
                return std::nullopt;
            }

            if (files.dir_status) {
                files.current_file.reset();
            } else {
                break;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief get the next sequence from the stream
     * @return sequence object, std::nullopt if end of stream
     */
    std::optional<sequence> next_sequence() {
        while (files.current_file || files.open_next_file()) {
            auto& fp = files.current_file;
            int l = fp->kseq_read_any();
            if (l == 0 && !preserve_empty_sequences_) {
                note_empty_sequence(*fp);
                continue;
            }
            if (l >= 0) {
                sequence r;
                r.id = fp->name_s();
                r.seq = fp->seq_s();
                r.is_fastq = fp->is_fastq_record();

                if (const char* c = fp->comment_s()){
                    r.comment = c;
                }
                if (const char* p = fp->plus_s()) {
                    r.plus = p;
                }
                if (r.is_fastq){
                    if (const char* q = fp->qual_s()) {
                        r.qual = q;
                    }
                }
                return r;
            }
            int gzip_code = Z_OK;
            std::string gzip_detail;
            if (l < 0 && fp->gzip_read_failed(gzip_code, gzip_detail)) {
                read_error_code_ = l < -1 ? l : -3;
                read_error_path_ = files.current_file_path();
                read_error_detail_ = "gzip/zlib code " +
                                     std::to_string(gzip_code);
                if (!gzip_detail.empty()) {
                    read_error_detail_ += ": " + gzip_detail;
                }
                return std::nullopt;
            }
            int pigz_status = 0;
            std::string pigz_detail;
            if (l < 0 && fp->pigz_read_failed(pigz_status, pigz_detail)) {
                read_error_code_ = -3;
                read_error_path_ = files.current_file_path();
                read_error_detail_ = pigz_detail;
                return std::nullopt;
            }
            if (l < -1) {
                read_error_code_ = l;
                read_error_path_ = files.current_file_path();
                return std::nullopt;
            }

            // Move to next file if in directory mode
            if (files.dir_status) {
                files.current_file.reset();
            } else {
                break;
            }
        }
        return std::nullopt;
    }
};

/**
 * @brief Memory utilities for read streaming, designed for monitoring memory usage of read streaming
 */
namespace readmem_utils {

    /**
     * @brief Get the memory bytes used by a string
     * @param s Input string
     * @return Number of bytes used by the string
     */
    static inline size_t bytes_of(const std::string& s) noexcept { 
        return s.capacity(); 
    }

    /**
     * @brief Get the memory bytes used by a read_streaming::sequence
     * @param r Input read_streaming::sequence
     * @return Number of bytes used by the sequence
     */
    static inline size_t bytes_of_rec(const read_streaming::sequence& r) noexcept {
        return bytes_of(r.id) + bytes_of(r.comment) + bytes_of(r.seq) +
               bytes_of(r.plus) + bytes_of(r.qual) + sizeof(r);
    }

    /**
     * @brief Get the memory bytes used by a vector of read_streaming::sequence
     * @tparam Vec Vector type
     * @param v Input vector
     * @return Number of bytes used by the vector
     */
    template <class Vec>
    static inline size_t bytes_of_vec(const Vec& v) noexcept {
        using T = typename Vec::value_type;
        size_t total = v.capacity() * sizeof(T);     // vector buffer itself
        for (auto const& e : v) total += bytes_of_rec(e); // payloads (strings' capacities)
        return total;
    }

    /**
     * @brief In-flight memory tracking for read streaming, tracks bytes and chunks in progress with `std::atomic<size_t>` counters
     * @param bytes Number of bytes in flight
     * @param chunks Number of chunks in flight
     * @param peak_bytes Peak bytes in flight
     * @param peak_chunks Peak chunks in flight
     */
    struct data_in_flight {
        std::atomic<size_t> bytes{0};
        std::atomic<size_t> chunks{0};
        std::atomic<size_t> peak_bytes{0};
        std::atomic<size_t> peak_chunks{0};

        void add(size_t b) noexcept {
            const size_t nb = bytes.fetch_add(b, std::memory_order_relaxed) + b;
            const size_t nc = chunks.fetch_add(1, std::memory_order_relaxed) + 1;
            size_t pb = peak_bytes.load(std::memory_order_relaxed);
            while (nb > pb && !peak_bytes.compare_exchange_weak(pb, nb, std::memory_order_relaxed)) {}
            size_t pc = peak_chunks.load(std::memory_order_relaxed);
            while (nc > pc && !peak_chunks.compare_exchange_weak(pc, nc, std::memory_order_relaxed)) {}
        }

        void sub(size_t b) noexcept {
            bytes.fetch_sub(b, std::memory_order_relaxed);
            chunks.fetch_sub(1, std::memory_order_relaxed);
        }
    };

    inline data_in_flight& inflight() { 
        static data_in_flight g; 
        return g; 
    }

    //
    struct scoped_in_flight {
        size_t b{0}; bool active{false};
        explicit scoped_in_flight(size_t bytes) : b(bytes), active(true) { inflight().add(b); }
        ~scoped_in_flight(){ if (active) inflight().sub(b); }
        scoped_in_flight(const scoped_in_flight&) = delete;
        scoped_in_flight& operator=(const scoped_in_flight&) = delete;
    };

    static inline double rss_gib() {
    #if defined(__APPLE__)
        // macOS: use Mach task info
        mach_task_basic_info info;
        mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) != KERN_SUCCESS) {
            return -1.0;
        }
        const double rss_bytes = static_cast<double>(info.resident_size);
        return rss_bytes / (1024.0 * 1024.0 * 1024.0);
    #elif defined(__linux__)
        // Linux: read statm
        FILE* f = std::fopen("/proc/self/statm", "r");
        if (!f) return -1.0;
        long long size_pages=0, res_pages=0;
        if (std::fscanf(f, "%lld %lld", &size_pages, &res_pages) != 2) { std::fclose(f); return -1.0; }
        std::fclose(f);
        const double page_kb = double(::sysconf(_SC_PAGESIZE)) / 1024.0;
        const double rss_kb = res_pages * page_kb;
        return rss_kb / (1024.0 * 1024.0);
    #else
        return -1.0;
    #endif
    }

    // --- tiny logger ---
    static inline void log_chunk(std::size_t chunk_bytes, std::size_t chunk_id=0) {
        const double mb = double(chunk_bytes) / (1024.0*1024.0);
        const double inflight_gib = double(inflight().bytes.load(std::memory_order_relaxed)) / (1024.0*1024.0*1024.0);
        const std::size_t inflight_chunks = inflight().chunks.load(std::memory_order_relaxed);
        const double rss = rss_gib();
        std::ostringstream message;
        message << "[radmem] chunk=" << chunk_id
                << " size_mb=" << std::fixed << std::setprecision(2) << mb
                << " inflight_chunks=" << inflight_chunks
                << " inflight_gib=" << std::setprecision(3) << inflight_gib
                << " rss_gib=" << rss << "\n";
        RAD_CERR << message.str();
    }
}


/**
 * @brief Chunk-based streaming with parallel decompression
 * 
 * Processes input files in chunks, using pigz for parallel decompression
 * and OpenMP tasks for parallel processing.
 */
template<typename T, typename Function> 
class chunk_streaming {
public:
    explicit chunk_streaming(size_t chunk_size, int pigz_threads = 4, int max_in_flight = -1) 
        : chunk_size_(chunk_size), pigz_threads_(pigz_threads), 
          max_in_flight_(max_in_flight), seqs_processed_(0), peak_in_flight_(0) {}

    void modify_chunk_size(size_t new_size) {
        chunk_size_ = new_size;
    }

    size_t get_chunk_size() const {
        return chunk_size_;
    }

    int get_mif() const {
        return max_in_flight_;
    }
  
    void process_chunks(const std::string& input_path, 
                    Function chunk_func, 
                    int num_threads, 
                    int64_t max_seqs = -1,
                    bool verbose = false,
                    size_t log_every_chunks = 0,
                    std::function<void()> boundary_poll = {}) {
        // Chunk acquisition and dispatch are intentionally serial. Parallel
        // work belongs inside the per-read processor; an outer OpenMP team
        // multiplied memory use and allowed exceptions to escape a structured
        // OpenMP block. This serial boundary is also the only safe place for
        // host callbacks such as R's user-interrupt check.
        (void)num_threads;
        file_streaming files(input_path, pigz_threads_);
        read_streaming reader(files);
        size_t seqs_processed = 0;
        size_t chunks_seen = 0;
        bool continue_processing = true;

        while (continue_processing) {
            if (boundary_poll) boundary_poll();

            std::vector<T> chunk;
            chunk.reserve(chunk_size_);
            while (chunk.size() < chunk_size_) {
                if (max_seqs >= 0 &&
                    seqs_processed >= static_cast<size_t>(max_seqs)) {
                    continue_processing = false;
                    break;
                }
                auto record = reader.next_sequence();
                if (!record) {
                    continue_processing = false;
                    break;
                }
                chunk.push_back(std::move(*record));
                ++seqs_processed;
            }

            if (chunk.empty()) break;

            const size_t chunk_bytes = readmem_utils::bytes_of_vec(chunk);
            const size_t chunk_id =
                next_chunk_id_.fetch_add(1, std::memory_order_relaxed);
            readmem_utils::scoped_in_flight guard(chunk_bytes);
            readmem_utils::log_chunk(chunk_bytes, chunk_id);

            if constexpr (std::is_same_v<
                              decltype(chunk_func(chunk, input_path)), bool>) {
                if (!chunk_func(chunk, input_path)) continue_processing = false;
            } else {
                chunk_func(chunk, input_path);
            }

            if (verbose) {
                ++chunks_seen;
                if (log_every_chunks == 0 ||
                    (chunks_seen % log_every_chunks == 0)) {
                    RAD_CERR << "[chunk_streaming] chunk=" << chunk_id
                             << " size=" << chunk_bytes
                             << " seqs=" << chunk.size() << "\n";
                }
            }
        }

        // kseq uses -1 for clean EOF and values below -1 for malformed or
        // truncated input. gzip/zlib integrity failures are normalized to -3.
        if (reader.read_failed()) {
            std::string message =
                "sequence read failed (kseq code " +
                std::to_string(reader.read_error_code()) + ") in " +
                reader.read_error_path();
            if (!reader.read_error_detail().empty()) {
                message += " (" + reader.read_error_detail() + ")";
            }
            throw std::runtime_error(message);
        }
        if (boundary_poll) boundary_poll();
        }


    private:
    size_t chunk_size_;
    int pigz_threads_;
    int max_in_flight_;  // -1 = unlimited, >0 = limited
    std::atomic<size_t> seqs_processed_;
    std::atomic<int> peak_in_flight_;
    std::atomic<size_t> next_chunk_id_{1};
};

/**
 * @brief pigz writing for parallel compression
 */
class pigz_writing {
public:
    inline FILE* open_pigz_pipe(const std::string& out_path, int threads, int level=1) {
        // Resolve pigz (PATH-first; honors RAD_PIGZ / RAD_NO_PIGZ). "" => skip.
        const std::string& pigz = pigz_reading::resolve_pigz();
        if (pigz.empty()) return nullptr;

        std::string cmd = pigz_reading::shell_quote_posix(pigz) +
                          " -c -p " + std::to_string(threads) + " -" +
                          std::to_string(level) + " > " +
                          pigz_reading::shell_quote_posix(out_path);
        FILE* fp = popen(cmd.c_str(), "w");
        if (!fp){
            return nullptr;
        } 
            // 16MB stdio buffer to reduce write() syscalls

            // removed because in debug mode, the buffer gets shared between all of the thread outputs
            // so then everything gets scrambled
            // TSPMO
        //static std::vector<char> buf(16u << 20);
        //setvbuf(fp, buf.data(), _IOFBF, buf.size());
        return fp;
    }

    inline void pigz_write(FILE* fp, const char* data, size_t n) {
        if (!n) return;
        size_t w = std::fwrite(data, 1, n, fp);
        if (w != n) {
            throw std::runtime_error("pigz_write: partial write");
        }
    }

    inline int close_pigz_pipe(FILE* fp) {
        if (!fp) return -1;
        std::fflush(fp);
        return pclose(fp);
    }
};

class sigstring_writing {
public:
    enum class format { 
        FASTQA, 
        SIGSTRING, 
        CSV 
    };

private:
    // Common
    format fmt;
    bool compress_;
    bool use_pigz_ = false;
    bool rc_umi_ = true; // reverse-complement UMIs on reverse reads (FASTQA only)
    bool finalized_ = false;
    std::string output_path_;

    // --- pigz path ---
    pigz_writing pigz_;
    FILE* pigz_fp_ = nullptr;
    std::string pigz_buf_;
    size_t pigz_flush_threshold_ = (8u << 20); // 8MB

    // --- in-process zlib/plain path ---
    std::ofstream file_;
    std::unique_ptr<ogzstream> gzip_file_;
    std::ostream* out_ = nullptr;

    inline void pigz_flush_if_big_() {
        if (pigz_buf_.size() >= pigz_flush_threshold_) {
            pigz_.pigz_write(pigz_fp_, pigz_buf_.data(), pigz_buf_.size());
            pigz_buf_.clear();
        }
    }

    template <typename T>
    inline void write_one_(const T& item) {
        switch (fmt) {
            case format::FASTQA:
                if (use_pigz_) { pigz_buf_.append(item.to_fastqa(rc_umi_)); pigz_flush_if_big_(); }
                else           { (*out_) << item.to_fastqa(rc_umi_); }
                break;
            case format::SIGSTRING:
                if (use_pigz_) { pigz_buf_.append(item.to_sigstring()); pigz_buf_.push_back('\n'); pigz_flush_if_big_(); }
                else           { (*out_) << item.to_sigstring() << "\n"; }
                break;
            case format::CSV:
                if (use_pigz_) { pigz_buf_.append(item.to_csv()); pigz_flush_if_big_(); }
                else           { (*out_) << item.to_csv(); }
                break;
        }
    }

public:
    /// @param output_path   base filename ("myout"), will add ".gz" if compress==true and it isn't present
    /// @param output_format which of the three formats to emit
    /// @param compress      if true, writes gzipped output (pigz if available; else zlib/boost)
    /// @param append        if true, opens in append mode; otherwise truncates
    /// @param pigz_threads  pigz worker threads (ignored if not using pigz)
    /// @param level         gzip compression level [1..9]
    sigstring_writing(std::string output_path, format output_format, bool compress, bool append, int pigz_threads, int level=1)
        : fmt(output_format), compress_(compress)
    {
        std::ios_base::openmode mode =
            std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc);

        if (compress_) {
            // ensure .gz suffix for either pigz or zlib
            if (output_path.size() < 3 || output_path.substr(output_path.size()-3) != ".gz")
                output_path += ".gz";

            output_path_ = output_path;

            // Try pigz first
            pigz_fp_ = pigz_.open_pigz_pipe(output_path, pigz_threads, level);
            if (pigz_fp_) {
                use_pigz_ = true;
                // nothing else to open; pigz writes the file
                return;
            }

            // Fallback: gzstream writes through zlib directly, which keeps the
            // embedded R build independent of compiled Boost libraries.
            if (append) {
                throw std::runtime_error(
                    "Appending compressed RAD output is not supported by the embedded backend");
            }
            gzip_file_ = std::make_unique<ogzstream>(output_path.c_str());
            if (!gzip_file_ || !*gzip_file_) {
                throw std::runtime_error("Failed to open file: " + output_path);
            }
            out_ = gzip_file_.get();
        } else {
            // Plain text (no compression)
            output_path_ = output_path;
            file_.open(output_path, mode);
            if (!file_.is_open()) throw std::runtime_error("Failed to open file: " + output_path);
            out_ = &file_;
        }
    }

    /// FASTQA only: if false, UMIs are written exactly as extracted (minus-strand
    /// on reverse reads) instead of being flipped to plus-strand. Mirrors the
    /// primary buffered path (SigString::to_fastqa_append) so --no-umi-rc is
    /// honored by the debug FASTQA output too.
    void set_rc_umi(bool v) { rc_umi_ = v; }

    template<typename T>
    void operator()(const std::vector<T>& chunk) {
        if (chunk.empty()) return;
        if (finalized_) {
            throw std::runtime_error("write attempted after finalizing: " +
                                     output_path_);
        }
        if (use_pigz_) {
            // Build big uncompressed buffer for pigz (few syscalls, fewer context switches)
            for (auto const& item : chunk) write_one_(item);
        } else {
            // Direct standard/gzip ostream path.
            for (auto const& item : chunk) write_one_(item);
            if (!out_ || !*out_) {
                throw std::runtime_error("output write failed: " +
                                         output_path_);
            }
        }
    }

    // Optional: adjust pigz flush size (in bytes)
    void set_pigz_flush_threshold(size_t bytes) { pigz_flush_threshold_ = bytes; }

    void finalize() {
        if (finalized_) return;
        finalized_ = true;
        std::exception_ptr failure;

        if (use_pigz_) {
            try {
                if (!pigz_buf_.empty()) {
                    pigz_.pigz_write(pigz_fp_, pigz_buf_.data(),
                                     pigz_buf_.size());
                    pigz_buf_.clear();
                }
            } catch (...) {
                failure = std::current_exception();
            }
            if (pigz_fp_) {
                if (std::fflush(pigz_fp_) != 0 && !failure) {
                    failure = std::make_exception_ptr(std::runtime_error(
                        "output flush failed: " + output_path_));
                }
                const int status = pclose(pigz_fp_);
                pigz_fp_ = nullptr;
                if (status != 0 && !failure) {
                    failure = std::make_exception_ptr(std::runtime_error(
                        "pigz close failed for " + output_path_ +
                        " (status " + std::to_string(status) + ")"));
                }
            }
            use_pigz_ = false;
        } else {
            if (out_) {
                out_->flush();
                if (!*out_) {
                    failure = std::make_exception_ptr(std::runtime_error(
                        "output flush failed: " + output_path_));
                }
            }
            if (gzip_file_) {
                gzip_file_->close();
                if (!*gzip_file_ && !failure) {
                    failure = std::make_exception_ptr(std::runtime_error(
                        "gzip close failed: " + output_path_));
                }
                gzip_file_.reset();
            }
            if (file_.is_open()) {
                file_.close();
                if (file_.fail() && !failure) {
                    failure = std::make_exception_ptr(std::runtime_error(
                        "output close failed: " + output_path_));
                }
            }
            out_ = nullptr;
        }

        if (failure) std::rethrow_exception(failure);
    }

    ~sigstring_writing() noexcept {
        try {
            finalize();
        } catch (...) {
        }
    }

    // in sigstring_writing (public)
    void operator()(const std::string& slab) {
        if (slab.empty()) return;
        if (finalized_) {
            throw std::runtime_error("write attempted after finalizing: " +
                                     output_path_);
        }
        if (use_pigz_) {
            pigz_buf_.append(slab);
            pigz_flush_if_big_();
        } else {
            out_->write(slab.data(), slab.size());
            if (!*out_) {
                throw std::runtime_error("output write failed: " +
                                         output_path_);
            }
        }
    }

     // reuses the existing std::string overload
    void operator()(const std::vector<std::string>& slabs) {
        if (slabs.empty()) return;
        for (auto const& s : slabs) {
            (*this)(s); 
        }
    }

};

/**
 * @brief Parallel writer for asynchronous output--works by queuing jobs and processing them in a separate thread.
 * Each write job is a function to perform the write operation, along with metadata.
 * This class uses a bounded standard queue and a dedicated writer thread.
 * 
 * What that means is that the main thread lines up write jobs and the writer thread executes them in the background
 * while the main thread continues processing. When single-threaded, the main thread does the writing itself. 
 * 
 * The standard queue keeps the embedded package independent of Boost.Lockfree,
 * which is not shipped by the BH R package.
 */
class parallel_writer {
public:
/**
 * @brief Constructs a parallel_writer and starts the writer thread
 */
    parallel_writer() : running(true), jobs_processed(0) {
        writer_thread = std::thread([this]() {
            this->process_jobs();
        });
    }
    
    ~parallel_writer() {
        try {
            stop();
        } catch (...) {
            // Destructors must not propagate.  The normal explicit stop() in
            // SigString::sigalign surfaces worker failures to the caller.
        }
    }
    
    /**
     * @brief Write multiple data items in parallel
     * @param fastqa_writer Reference to `sigstring_writing` object for fastqa output
     * @param data `vector` of data items to write
     */
    template<typename T>
    void write(sigstring_writing& fastqa_writer, std::vector<T>& data) {
        if (data.empty()) return;
        
        auto job_function = [
            fastq_ptr=&fastqa_writer,
            data_vec=std::move(data)
        ]() mutable {
            (*fastq_ptr)(data_vec);
            std::vector<T>().swap(data_vec);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data.size();
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        enqueue(job);
    }
    
    /**
     * @brief Write multiple string slabs in parallel
     * @param fastqa_writer Reference to sigstring_writing object for FASTQ/A output
     * @param slabs Vector of string slabs to write
     */
    void write_slabs(sigstring_writing& fastqa_writer, std::vector<std::string>&& slabs) {
        size_t total_bytes = 0;
        for (auto const& s : slabs) total_bytes += s.size();
        if (total_bytes == 0) return;
        
        auto* job = new WriteJob();
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = total_bytes;
        job->queued_time = std::chrono::high_resolution_clock::now();
        job->func = [fastq_ptr=&fastqa_writer, slabs_vec=std::move(slabs)]() mutable {
            (*fastq_ptr)(slabs_vec);
            std::vector<std::string>().swap(slabs_vec);
        };
        
        enqueue(job);
    }
    
    /**
     * @brief Write a raw string in parallel
     * @param writer Reference to `sigstring_writing` object for output
     * @param data `string` data to write
     */
    void write_raw_string(sigstring_writing& writer, std::string&& data) {
        if (data.empty()) return;
        
        size_t data_size = data.size();
        
        auto job_function = [
            writer_ptr=&writer,
            data_str=std::move(data)
        ]() mutable {
            (*writer_ptr)(data_str);
            std::string().swap(data_str);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data_size;
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        enqueue(job);
    }

    /**
     * @brief Write debug information in parallel
     * @param sig_writer Pointer to sigstring_writing for SIGSTRING output
     * @param csv_writer Pointer to sigstring_writing for CSV output
     * @param fastqa_writer Pointer to sigstring_writing for FASTQ/A output
     * @param data Vector of data to write
     */
    template<typename T>
    void write_debug(sigstring_writing* sig_writer, sigstring_writing* csv_writer, sigstring_writing* fastqa_writer, std::vector<T>& data) {
        if (data.empty()) return;
        
        auto job_function = [
            sig_ptr=sig_writer,
            csv_ptr=csv_writer,
            fastqa_ptr=fastqa_writer,
            data_vec=std::move(data)
        ]() mutable {
            (*sig_ptr)(data_vec);
            (*csv_ptr)(data_vec);
            (*fastqa_ptr)(data_vec);
            std::vector<T>().swap(data_vec);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data.size();
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        enqueue(job);
    }

    void stop() {
        if (stop_completed) return;
        running.store(false, std::memory_order_release);
        queue_not_empty.notify_all();
        queue_not_full.notify_all();
        
        if (writer_thread.joinable()) {
            writer_thread.join();
        }
        stop_completed = true;
        std::exception_ptr error = writer_error;
        writer_error = nullptr;
        
        RAD_COUT << "\n[Writer Thread] Statistics:" << std::endl;
        RAD_COUT << "  Total jobs processed: " << jobs_processed.load() << std::endl;
        RAD_COUT << "  Total write time: " << total_write_time_ms.load() / 1000.0 << " seconds" << std::endl;
        if (jobs_processed.load() > 0) {
            RAD_COUT << "  Average write time per chunk: " << total_write_time_ms.load() / jobs_processed.load() << " ms" << std::endl;
            RAD_COUT << "  Average queue time per chunk: " << total_queue_time_ms.load() / jobs_processed.load() << " ms" << std::endl;
        }
        if (error) std::rethrow_exception(error);
    }

    double get_total_write_time_ms() const {
        return total_write_time_ms.load();
    }
    
    size_t get_jobs_processed() const {
        return jobs_processed.load();
    }

private:
/**
 * @brief Structure representing a write job
 * @param func Function to perform the write operation
 * @param chunk_id Unique chunk identifier for logging
 * @param size Size of data to write
 * @param queued_time Time when job was queued
 */
    struct WriteJob {
        std::function<void()> func;
        size_t chunk_id;
        size_t size;
        std::chrono::time_point<std::chrono::high_resolution_clock> queued_time;
    };
    
    static constexpr size_t max_queued_jobs = 64;
    std::queue<WriteJob*> write_queue;
    std::mutex queue_mutex;
    std::condition_variable queue_not_empty;
    std::condition_variable queue_not_full;
    std::exception_ptr writer_error;
    bool stop_completed = false;
    std::atomic<bool> running; /// we writing? y/n
    std::atomic<size_t> next_chunk_id{1}; /// next chunk ID to assign
    std::atomic<size_t> jobs_processed{0}; /// total jobs processed
    std::thread writer_thread; /// writer thread
    
    std::atomic<double> total_write_time_ms{0.0};
    std::atomic<double> total_queue_time_ms{0.0};
    
    /**
     * @brief writer function per thread to process queued jobs
     */
    void enqueue(WriteJob* job) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        queue_not_full.wait(lock, [this]() {
            return write_queue.size() < max_queued_jobs ||
                   !running.load(std::memory_order_acquire);
        });
        if (!running.load(std::memory_order_acquire)) {
            delete job;
            if (writer_error) std::rethrow_exception(writer_error);
            throw std::runtime_error("RAD output writer stopped unexpectedly");
        }
        write_queue.push(job);
        lock.unlock();
        queue_not_empty.notify_one();
    }

    void process_jobs() {
        while (true) {
            WriteJob* job = nullptr;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                queue_not_empty.wait(lock, [this]() {
                    return !write_queue.empty() ||
                           !running.load(std::memory_order_acquire);
                });
                if (write_queue.empty() &&
                    !running.load(std::memory_order_acquire)) {
                    break;
                }
                job = write_queue.front();
                write_queue.pop();
            }
            queue_not_full.notify_one();

            try {
                auto start_time = std::chrono::high_resolution_clock::now();
                auto queue_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    start_time - job->queued_time).count();
                total_queue_time_ms.store(
                    total_queue_time_ms.load(std::memory_order_relaxed) + queue_time_ms,
                    std::memory_order_relaxed);

                job->func();

                auto end_time = std::chrono::high_resolution_clock::now();
                auto write_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    end_time - start_time).count();
                total_write_time_ms.store(
                    total_write_time_ms.load(std::memory_order_relaxed) + write_time_ms,
                    std::memory_order_relaxed);
                
                jobs_processed.fetch_add(1, std::memory_order_relaxed);
                
                const double size_mib = double(job->size) / (1024.0 * 1024.0);
                RAD_COUT << "\n[io_writer] chunk=" << job->chunk_id
                          << " size_mib=" << std::fixed << std::setprecision(2) << size_mib
                          << " queued_ms=" << queue_time_ms
                          << " wrote_ms=" << write_time_ms
                          << std::endl;
                
                delete job;
            } catch (...) {
                delete job;
                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    writer_error = std::current_exception();
                    running.store(false, std::memory_order_release);
                }
                queue_not_empty.notify_all();
                queue_not_full.notify_all();
                break;
            }
        }

        // Drop only jobs that could not run after a writer error.  A normal
        // stop reaches the loop only after all queued work is drained.
        if (writer_error) {
            std::lock_guard<std::mutex> lock(queue_mutex);
            while (!write_queue.empty()) {
                delete write_queue.front();
                write_queue.pop();
            }
        }
    }
};

#pragma once

#ifdef RRAD_EMBEDDED

#include <cstdarg>
#include <cstdio>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

// The embedded RAD engine cannot write through stdout/stderr or the standard
// C++ console streams: R may replace those connections, and R CMD check treats
// direct use from package code as non-portable.  Keep output in a package-owned
// sink instead.  The sink is process-local, supports nested captures, and
// serializes writes from RAD's worker threads without calling the R API there.
namespace rrad_embedded_console {

struct CaptureBuffer {
    std::string text;
};

inline std::recursive_mutex& capture_mutex() {
    static std::recursive_mutex value;
    return value;
}

inline std::vector<CaptureBuffer*>& capture_stack() {
    static std::vector<CaptureBuffer*> value;
    return value;
}

class ConsoleStream {
public:
    template <typename T>
    ConsoleStream& operator<<(const T& value) {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        formatter_ << value;
        publish_formatted_text();
        return *this;
    }

    using OstreamManipulator = std::ostream& (*)(std::ostream&);
    using IosManipulator = std::ios& (*)(std::ios&);
    using IosBaseManipulator = std::ios_base& (*)(std::ios_base&);

    ConsoleStream& operator<<(OstreamManipulator manipulator) {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        manipulator(formatter_);
        publish_formatted_text();
        return *this;
    }

    ConsoleStream& operator<<(IosManipulator manipulator) {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        manipulator(formatter_);
        publish_formatted_text();
        return *this;
    }

    ConsoleStream& operator<<(IosBaseManipulator manipulator) {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        manipulator(formatter_);
        publish_formatted_text();
        return *this;
    }

private:
    void publish_formatted_text() {
        const std::string addition = formatter_.str();
        formatter_.str("");
        formatter_.clear();
        if (!addition.empty() && !capture_stack().empty()) {
            capture_stack().back()->text += addition;
        }
    }

    std::ostringstream formatter_;
};

inline ConsoleStream& out() {
    static ConsoleStream value;
    return value;
}

inline ConsoleStream& err() {
    static ConsoleStream value;
    return value;
}

class ScopedCapture {
public:
    ScopedCapture() {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        capture_stack().push_back(&buffer_);
    }

    ScopedCapture(const ScopedCapture&) = delete;
    ScopedCapture& operator=(const ScopedCapture&) = delete;

    ~ScopedCapture() noexcept {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        auto& stack = capture_stack();
        for (auto it = stack.end(); it != stack.begin();) {
            --it;
            if (*it == &buffer_) {
                stack.erase(it);
                break;
            }
        }
    }

    std::string str() const {
        std::lock_guard<std::recursive_mutex> guard(capture_mutex());
        return buffer_.text;
    }

private:
    CaptureBuffer buffer_;
};

inline void errorf(const char* format, ...) {
    if (!format) return;

    va_list arguments;
    va_start(arguments, format);
    va_list size_arguments;
    va_copy(size_arguments, arguments);
    const int required = std::vsnprintf(nullptr, 0, format, size_arguments);
    va_end(size_arguments);

    if (required > 0) {
        std::vector<char> message(static_cast<std::size_t>(required) + 1U);
        (void)std::vsnprintf(message.data(), message.size(), format, arguments);
        err() << message.data();
    }
    va_end(arguments);
}

}  // namespace rrad_embedded_console

#define RAD_COUT (::rrad_embedded_console::out())
#define RAD_CERR (::rrad_embedded_console::err())

#else

#include <iostream>

#define RAD_COUT std::cout
#define RAD_CERR std::cerr

#endif

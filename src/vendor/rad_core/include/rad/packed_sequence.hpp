#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

__extension__ typedef unsigned __int128 uint128_t;

template <typename Word>
constexpr Word packed_mask(std::size_t base_length) {
    static_assert(std::is_same_v<Word, uint64_t> || std::is_same_v<Word, uint128_t>,
                  "Packed DNA requires an unsigned 64-bit or 128-bit word");
    constexpr std::size_t capacity = sizeof(Word) * 4;
    if (base_length > capacity)
        throw std::length_error("Packed DNA length exceeds the selected word width");
    if (base_length == 0) return Word{0};
    if (base_length == capacity) return ~Word{0};
    return (Word{1} << (2 * base_length)) - Word{1};
}

namespace packed_sequence_detail {
inline uint64_t mix_word(uint64_t value) noexcept {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

inline std::size_t hash_words(uint64_t low, uint64_t high, uint16_t length,
                              bool initialized = true) noexcept {
    uint64_t value = mix_word(low);
    value ^= mix_word(high + UINT64_C(0x243f6a8885a308d3));
    value ^= mix_word(static_cast<uint64_t>(length) |
                      (static_cast<uint64_t>(initialized) << 16));
    if constexpr (sizeof(std::size_t) < sizeof(uint64_t))
        value ^= value >> 32;
    return static_cast<std::size_t>(value);
}

template <typename Word>
inline std::size_t hash_word(Word value, uint16_t length) noexcept {
    uint64_t high = 0;
    if constexpr (sizeof(Word) > sizeof(uint64_t))
        high = static_cast<uint64_t>(value >> 64);
    return hash_words(static_cast<uint64_t>(value), high, length);
}

inline uint8_t nucleotide_code(char base) noexcept {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'T': return 2;
        case 'G': return 3;
        default: return 4;
    }
}

template <typename Word>
bool encode_word(std::string_view sequence, Word& value) noexcept {
    value = 0;
    for (char base : sequence) {
        const uint8_t code = nucleotide_code(base);
        if (code > 3) return false;
        value = (value << 2) | static_cast<Word>(code);
    }
    return true;
}

template <typename Word>
std::string decimal_word(Word value) {
    if (value == 0) return "0";
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + static_cast<unsigned>(value % 10)));
        value /= 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}
}

template <typename Word>
struct packed_barcode {
    static_assert(std::is_same_v<Word, uint64_t> || std::is_same_v<Word, uint128_t>,
                  "Packed barcodes require an unsigned 64-bit or 128-bit word");
    using word_type = Word;
    Word bits = 0;
    uint16_t length = 0;

    constexpr bool operator==(const packed_barcode& other) const noexcept {
        return length == other.length && bits == other.bits;
    }
    constexpr bool operator!=(const packed_barcode& other) const noexcept {
        return !(*this == other);
    }
    constexpr bool is_valid() const noexcept {
        return length > 0 && length <= sizeof(Word) * 4;
    }
};

class int64_seq {
    union sequence_word {
        uint64_t narrow;
        uint128_t wide;
        constexpr sequence_word() : narrow(0) {}
    };
    sequence_word word_;

public:
    uint16_t length;

private:
    bool initialized_;

    void set_invalid() noexcept {
        word_.narrow = 0;
        length = 0;
        initialized_ = true;
    }

    template <typename Word>
    void assign_word(Word raw, uint16_t seq_length) {
        static_assert(std::is_integral_v<Word> || std::is_same_v<Word, uint128_t>,
                      "Packed DNA requires an integer word");
        static_assert(sizeof(Word) <= sizeof(uint128_t), "Packed DNA word exceeds 128 bits");
        if (seq_length > 64)
            throw std::length_error("RAD packed DNA supports sequences of at most 64 bases");
        if (seq_length <= 32) {
            if constexpr (sizeof(Word) > sizeof(uint64_t)) {
                if ((static_cast<uint128_t>(raw) >> 64) != 0)
                    throw std::overflow_error("Packed DNA word exceeds its 64-bit sequence storage");
            }
            word_.narrow = static_cast<uint64_t>(raw);
        } else {
            if constexpr (sizeof(Word) <= sizeof(uint64_t))
                word_.wide = static_cast<uint128_t>(static_cast<uint64_t>(raw));
            else
                word_.wide = static_cast<uint128_t>(raw);
        }
        length = seq_length;
        initialized_ = true;
    }

public:
    int64_seq() noexcept : word_{}, length(0), initialized_(false) {}

    explicit int64_seq(uint64_t raw, uint16_t seq_length) : int64_seq() {
        assign_word(raw, seq_length);
    }

    explicit int64_seq(std::string_view sequence) : int64_seq() {
        sequence_to_bits(sequence);
    }

    template <typename Word>
    explicit int64_seq(const rrad_rad_core::packed_barcode<Word>& sequence) : int64_seq() {
        assign_word(sequence.bits, sequence.length);
    }

    template <typename Word>
    static int64_seq from_word(Word raw, uint16_t seq_length) {
        int64_seq result;
        result.assign_word(raw, seq_length);
        return result;
    }

    void sequence_to_bits(std::string_view sequence) {
        if (sequence.size() > 64)
            throw std::length_error("RAD packed DNA supports sequences of at most 64 bases");
        const uint16_t seq_length = static_cast<uint16_t>(sequence.size());
        if (seq_length <= 32) {
            uint64_t raw = 0;
            if (!packed_sequence_detail::encode_word(sequence, raw)) {
                set_invalid();
                return;
            }
            word_.narrow = raw;
        } else {
            uint128_t raw = 0;
            if (!packed_sequence_detail::encode_word(sequence, raw)) {
                set_invalid();
                return;
            }
            word_.wide = raw;
        }
        length = seq_length;
        initialized_ = true;
    }

    bool empty() const noexcept { return !initialized_; }
    bool is_valid() const noexcept { return initialized_ && length > 0 && length <= 64; }
    bool is_wide() const noexcept { return length > 32; }

    uint64_t as_uint64() const {
        if (is_wide())
            throw std::length_error("A sequence longer than 32 bases requires 128-bit storage");
        return word_.narrow;
    }

    uint128_t as_uint128() const noexcept {
        return is_wide() ? word_.wide : static_cast<uint128_t>(word_.narrow);
    }

    template <typename Consumer>
    decltype(auto) visit(Consumer&& consume) const {
        if (is_wide()) return std::forward<Consumer>(consume)(word_.wide);
        return std::forward<Consumer>(consume)(word_.narrow);
    }

    uint8_t base_at(std::size_t index) const {
        if (index >= length) throw std::out_of_range("Packed DNA base index is out of range");
        const std::size_t shift = 2 * (length - 1 - index);
        if (is_wide()) return static_cast<uint8_t>((word_.wide >> shift) & uint128_t{3});
        return static_cast<uint8_t>((word_.narrow >> shift) & uint64_t{3});
    }

    int64_seq slice(std::size_t start, std::size_t count) const {
        if (start > length || count > length - start)
            throw std::out_of_range("Packed DNA slice is out of range");
        if (count == 0)
            return empty() ? int64_seq{} : int64_seq(uint64_t{0}, 0);
        const std::size_t shift = 2 * (length - start - count);
        return visit([&](auto raw) {
            using Word = decltype(raw);
            return from_word(static_cast<Word>((raw >> shift) & packed_mask<Word>(count)),
                             static_cast<uint16_t>(count));
        });
    }

    std::string bits_to_sequence() const {
        static constexpr char bases[] = {'A', 'C', 'T', 'G'};
        std::string sequence(length, 'A');
        for (std::size_t i = 0; i < length; ++i)
            sequence[i] = bases[base_at(i)];
        return sequence;
    }

    bool operator==(const int64_seq& other) const noexcept {
        if (length != other.length || initialized_ != other.initialized_) return false;
        return is_wide() ? word_.wide == other.word_.wide : word_.narrow == other.word_.narrow;
    }
    bool operator!=(const int64_seq& other) const noexcept { return !(*this == other); }

    std::size_t hash_value() const noexcept {
        if (empty()) return packed_sequence_detail::hash_words(0, 0, length, false);
        if (is_wide()) return packed_sequence_detail::hash_word(word_.wide, length);
        return packed_sequence_detail::hash_word(word_.narrow, length);
    }

    std::size_t storage_bytes() const noexcept { return sizeof(*this); }

    void print_int64_seq() const {
        RAD_COUT << "Sequence length: " << length << "\nPacked bits: ";
        visit([](auto raw) { RAD_COUT << packed_sequence_detail::decimal_word(raw); });
        RAD_COUT << '\n';
    }
};

inline bool parse_packed_decimal(std::string_view token, uint16_t length,
                                 int64_seq& out, bool strict = false) {
    if (token.empty() || length > 64 || (strict && length == 0)) return false;
    size_t position = 0;
    while (position < token.size() &&
           std::isspace(static_cast<unsigned char>(token[position]))) ++position;
    if (position == token.size()) return false;
    if (token[position] == '-') {
        try {
            size_t consumed = 0;
            const uint64_t raw = static_cast<uint64_t>(std::stoll(std::string(token), &consumed));
            if (strict && (consumed != token.size() ||
                (length < 32 && (raw & ~packed_mask<uint64_t>(length)) != 0))) return false;
            out = int64_seq::from_word(raw, length);
            return true;
        } catch (const std::exception&) { return false; }
    }
    if (token[position] == '+') ++position;
    if (position == token.size() || token[position] < '0' || token[position] > '9') return false;
    auto parse_unsigned = [&](auto zero) {
        using Word = decltype(zero);
        Word raw = 0;
        const Word maximum = ~Word{0};
        size_t end = position;
        for (; end < token.size() && token[end] >= '0' && token[end] <= '9'; ++end) {
            const Word digit = static_cast<Word>(token[end] - '0');
            if (raw > (maximum - digit) / 10) return false;
            raw = raw * 10 + digit;
        }
        if (strict && (end != token.size() || (raw & ~packed_mask<Word>(length)) != 0)) return false;
        out = int64_seq::from_word(raw, length);
        return true;
    };
    if (length <= 32) return parse_unsigned(uint64_t{0});
    return parse_unsigned(uint128_t{0});
}

} // namespace rrad_rad_core
namespace std {
template <>
struct hash<rrad_rad_core::int64_seq> {
    std::size_t operator()(const rrad_rad_core::int64_seq& sequence) const noexcept {
        return sequence.hash_value();
    }
};
template <typename Word>
struct hash<rrad_rad_core::packed_barcode<Word>> {
    std::size_t operator()(const rrad_rad_core::packed_barcode<Word>& sequence) const noexcept {
        return rrad_rad_core::packed_sequence_detail::hash_word(sequence.bits, sequence.length);
    }
};
}

namespace rrad_rad_core {

#pragma once

namespace rad_discovery {
namespace layout_detail {

struct Part { std::string id, seq, klass; char poly = 0; };
struct Hit { int part, start, end, edits; bool reverse; };
struct Vocabulary {
    std::vector<Part> parts;
    std::map<std::pair<int, int>, std::vector<std::vector<Row>>> spans;
};

std::string csv_field(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) { if (c == '"') out += '"'; out += c; }
    return out + '"';
}

int median(std::vector<int> values) {
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

bool canonical(const std::string& seq) {
    return !seq.empty() && seq.find_first_not_of("ACGT") == std::string::npos;
}

int add_part(Vocabulary& v, const Row& row) {
    char poly = row.klass == "poly_t" ? 'T' : row.klass == "poly_a" ? 'A' : 0;
    std::string seq = poly ? std::string(12, poly) : row.sequence;
    if (!poly && (!canonical(seq) || seq.size() < 8)) return -1;
    for (int i = 0; i < static_cast<int>(v.parts.size()); ++i)
        if (v.parts[i].seq == seq) return i;
    v.parts.push_back({row.id, seq, row.klass, poly});
    return static_cast<int>(v.parts.size() - 1);
}

void load_source(Vocabulary& v, const std::string& path, bool parts_only) {
    csv::CSVFormat fmt;
    fmt.delimiter(',').header_row(ReadLayout::detect_header_row(path));
    csv::CSVReader reader(path, fmt);
    const auto headers = reader.get_col_names();
    auto field = [&](csv::CSVRow& r, const std::string& name) {
        if (std::find(headers.begin(), headers.end(), name) == headers.end()) return std::string();
        return seq_utils::trim(r[name].get<std::string>());
    };
    if (std::find(headers.begin(), headers.end(), "id") == headers.end() ||
        std::find(headers.begin(), headers.end(), "seq") == headers.end())
        throw std::runtime_error("Discovery parts require id and seq columns: " + path);
    std::vector<Row> variables;
    int previous = -1;
    for (auto& r : reader) {
        Row row;
        row.id = field(r, "id"); row.sequence = field(r, "seq");
        std::transform(row.sequence.begin(), row.sequence.end(), row.sequence.begin(), ::toupper);
        row.type = parts_only ? "static" : field(r, "type");
        row.klass = field(r, "class");
        row.whitelist = field(r, "whitelist"); row.flags = field(r, "flags");
        row.length_spec = field(r, "expected_length");
        auto lengths = ReadLayout::parse_length_candidates(row.length_spec);
        row.length = lengths.empty() ? 0 : lengths.back();
        if (row.id == "poly_t" || row.id == "poly_a") row.klass = row.id;
        const std::string direction = field(r, "direction");
        if (direction.rfind("reverse", 0) == 0) continue;
        if (row.type == "static") {
            int part = add_part(v, row);
            if (part < 0) continue;
            if (!parts_only) v.spans[{previous, part}].push_back(variables);
            previous = part;
            variables.clear();
        } else if (!parts_only) {
            variables.push_back(row);
        }
    }
    if (!parts_only) v.spans[{previous, -2}].push_back(variables);
}

Vocabulary vocabulary(const Options& options, std::vector<std::string>& notes) {
    Vocabulary v;
    config_utils::load_layout_overrides_once();
    std::vector<std::string> names;
    for (const auto& item : config_utils::layout_files) names.push_back(item.first);
    std::sort(names.begin(), names.end());
    std::set<std::string> loaded;
    for (const auto& name : names) {
        if (options.boundary_poll) options.boundary_poll();
        try {
            auto path = config_utils::get_read_layout(name);
            if (loaded.insert(path).second) load_source(v, path, false);
        } catch (const std::exception& e) {
            notes.push_back("Skipped layout vocabulary " + name + ": " + e.what());
        }
    }
    if (!options.parts_csv.empty()) load_source(v, options.parts_csv, true);
    if (v.parts.empty()) throw std::runtime_error("No usable static parts are available");
    return v;
}

std::vector<Hit> find_hits(const std::string& sequence, const Vocabulary& v, double error_ratio) {
    std::vector<Hit> hits;
    for (int p = 0; p < static_cast<int>(v.parts.size()); ++p) {
        const auto& part = v.parts[p];
        for (bool reverse : {false, true}) {
            std::string query = reverse ? seq_utils::revcomp(part.seq) : part.seq;
            if (part.poly) {
                for (int pos = 0; pos < static_cast<int>(sequence.size());) {
                    int end = pos;
                    while (end < static_cast<int>(sequence.size()) && sequence[end] == query[0]) ++end;
                    if (end - pos >= 12) hits.push_back({p, pos, end, 0, reverse});
                    pos = std::max(pos + 1, end);
                    if (hits.size() >= 256) break;
                }
                continue;
            }
            int budget = static_cast<int>(std::floor(query.size() * error_ratio));
            auto aligned = edlibAlign(query.data(), static_cast<int>(query.size()),
                                      sequence.data(), static_cast<int>(sequence.size()),
                                      edlibNewAlignConfig(budget, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0));
            if (aligned.status == EDLIB_STATUS_OK && aligned.editDistance >= 0) {
                for (int k = 0; k < std::min(aligned.numLocations, 32); ++k)
                    hits.push_back({p, aligned.startLocations[k], aligned.endLocations[k] + 1,
                                    aligned.editDistance, reverse});
            }
            edlibFreeAlignResult(aligned);
        }
    }
    // Alternative primer names and nested primer sequences must not become two parts at one site.
    std::sort(hits.begin(), hits.end(), [&](const Hit& a, const Hit& b) {
        int sa = static_cast<int>(v.parts[a.part].seq.size()) - 3 * a.edits;
        int sb = static_cast<int>(v.parts[b.part].seq.size()) - 3 * b.edits;
        if (sa != sb) return sa > sb;
        return std::tie(a.start, a.part, a.reverse) < std::tie(b.start, b.part, b.reverse);
    });
    std::vector<Hit> kept;
    for (const auto& hit : hits) {
        bool overlaps = false;
        for (const auto& other : kept) {
            if (other.reverse != hit.reverse) continue;
            int overlap = std::min(hit.end, other.end) - std::max(hit.start, other.start);
            if (overlap > 0 && overlap * 2 >= std::min(hit.end - hit.start, other.end - other.start)) {
                overlaps = true; break;
            }
        }
        if (!overlaps) kept.push_back(hit);
    }
    return kept;
}

struct Gap {
    int middle = 0;
    bool stable = false;
    std::size_t count = 0;
};

Gap describe_gap(const std::vector<int>& values) {
    Gap g;
    if (values.empty()) return g;
    g.middle = median(values); g.count = values.size();
    std::size_t near = 0;
    for (int value : values) near += std::abs(value - g.middle) <= 2;
    g.stable = near * 10 >= values.size() * 8;
    return g;
}

std::vector<Row> infer_gap(const Vocabulary& v, int left, int right, const Gap& gap,
                           LayoutInference& result, int& serial) {
    if (gap.count == 0) {
        if (left >= 0 && right >= 0) {
            result.complete = false;
            result.notes.push_back("No same-read spacing evidence for an adjacent pair of selected parts");
        }
        return {};
    }
    if (gap.middle <= 2) return {};
    auto known = v.spans.find({left, right});
    std::vector<Row> chosen;
    int best_distance = INT_MAX;
    bool role_tie = false;
    auto signature = [](const std::vector<Row>& rows) {
        std::string s;
        for (const auto& row : rows) s += row.klass + ":" + std::to_string(row.length) + ";";
        return s;
    };
    if (known != v.spans.end()) {
        for (const auto& span : known->second) {
            if (span.empty()) continue;
            int total = 0; bool insert = false;
            for (const auto& row : span) { total += row.length; insert |= row.klass == "read"; }
            int distance = insert ? (gap.middle >= total + 20 ? 0 : INT_MAX) : std::abs(total - gap.middle);
            if (!insert && distance > 6) continue;
            if (distance < best_distance) { chosen = span; best_distance = distance; role_tie = false; }
            else if (distance == best_distance && signature(chosen) != signature(span)) role_tie = true;
        }
    }
    if (!chosen.empty()) {
        if (role_tie) {
            result.complete = false;
            result.notes.push_back("Alternative variable-region contracts fit the same observed gap");
        }
        int total = 0, umi = -1;
        bool insert = false;
        for (int i = 0; i < static_cast<int>(chosen.size()); ++i) {
            total += chosen[i].length;
            if (chosen[i].klass == "umi") umi = i;
            insert |= chosen[i].klass == "read";
            if (chosen[i].length_spec.find('-') != std::string::npos) {
                result.complete = false;
                result.notes.push_back("Individual component lengths remain unresolved within a bounded gap");
            }
            // Catalog identity is detected from raw barcode observations, never inherited from a template.
            chosen[i].whitelist.clear();
        }
        if (!insert && gap.stable && umi >= 0) {
            int residual = chosen[umi].length + gap.middle - total;
            if (residual >= 4 && residual <= 32) {
                chosen[umi].length = residual;
                chosen[umi].length_spec = std::to_string(residual);
                result.notes.push_back("UMI role uses known flanking-part semantics; length is the observed gap residual");
            } else result.complete = false;
        } else if (!insert && (!gap.stable || total != gap.middle)) {
            result.complete = false;
            result.notes.push_back("Variable-region lengths are not resolved by a stable bounded gap");
        }
        if ((left == -1 || right == -2) && !insert) {
            result.complete = false;
            result.notes.push_back("Read-end distance does not establish an unbounded UMI or barcode endpoint");
        }
        return chosen;
    }
    Row row;
    row.id = "region_" + std::to_string(++serial); row.type = "variable";
    if (gap.middle >= 50 || !gap.stable) {
        row.klass = "read";
    } else {
        row.klass = "spacer"; row.length = gap.middle;
        row.length_spec = std::to_string(row.length);
        result.complete = false;
        result.notes.push_back("Fixed gap " + std::to_string(gap.middle) + " bp has no established barcode/UMI role");
    }
    return {row};
}

void make_recipes(LayoutInference& result) {
    for (std::size_t b = 0; b < result.rows.size(); ++b) {
        const auto& row = result.rows[b];
        if (row.klass != "barcode") continue;
        auto lengths = ReadLayout::parse_length_candidates(row.length_spec);
        if (lengths.empty() && row.length > 0) lengths.push_back(row.length);
        int left = static_cast<int>(b) - 1, right = static_cast<int>(b) + 1;
        std::vector<int> before{0}, after{0};
        auto extend = [](std::vector<int>& offsets, const Row& neighbor) {
            auto choices = ReadLayout::parse_length_candidates(neighbor.length_spec);
            if (choices.empty()) choices.push_back(neighbor.length);
            std::set<int> expanded;
            for (int offset : offsets) for (int length : choices)
                if (length > 0 && offset + length <= 256) expanded.insert(offset + length);
            offsets.assign(expanded.begin(), expanded.end());
        };
        while (left >= 0 && result.rows[left].type != "static") {
            if (result.rows[left].length <= 0) { left = -1; break; }
            extend(before, result.rows[left--]);
        }
        while (right < static_cast<int>(result.rows.size()) && result.rows[right].type != "static") {
            if (result.rows[right].length <= 0) { right = static_cast<int>(result.rows.size()); break; }
            extend(after, result.rows[right++]);
        }
        auto literal = [](const Row& r) {
            if (r.klass == "poly_t") return std::string(12, 'T');
            if (r.klass == "poly_a") return std::string(12, 'A');
            return r.sequence;
        };
        // Prefer a literal adapter over a poly run, whose edge can absorb barcode/UMI bases.
        bool use_left = left >= 0 && (right >= static_cast<int>(result.rows.size()) ||
            result.rows[left].klass.rfind("poly_", 0) != 0 || result.rows[right].klass.rfind("poly_", 0) == 0);
        for (int length : lengths) {
            const auto& offsets = use_left ? before : after;
            for (int offset : offsets) {
                Recipe recipe;
                recipe.id = row.id; recipe.length = length;
                if (use_left) recipe.pattern = literal(result.rows[left]) + std::string(offset, 'n') + std::string(length, 'N');
                else if (right < static_cast<int>(result.rows.size()))
                    recipe.pattern = std::string(length, 'N') + std::string(offset, 'n') + literal(result.rows[right]);
                if (!recipe.pattern.empty()) result.recipes.push_back(std::move(recipe));
            }
        }
    }
}

// Whitelist-only discovery remains useful when a read contains the barcode-side
// adapter but lacks the other boundaries needed to generate a complete layout.
void fallback_recipes(const Vocabulary& vocab, const std::vector<std::size_t>& support,
                      LayoutInference& result) {
    if (!result.recipes.empty()) return;
    std::set<std::tuple<std::string, int, int, int>> seen;
    int serial = 0;
    for (const auto& entry : vocab.spans) {
        int left = entry.first.first, right = entry.first.second;
        bool left_seen = left >= 0 && support[left] >= 3;
        bool right_seen = right >= 0 && support[right] >= 3;
        if (!left_seen && !right_seen) continue;
        for (const auto& span : entry.second) {
            LayoutInference candidate;
            auto append_part = [&](int part) {
                Row row;
                row.id = vocab.parts[part].id; row.sequence = vocab.parts[part].seq;
                row.type = "static"; row.klass = vocab.parts[part].poly ?
                    (vocab.parts[part].poly == 'T' ? "poly_t" : "poly_a") : vocab.parts[part].klass;
                candidate.rows.push_back(std::move(row));
            };
            if (left_seen) append_part(left);
            candidate.rows.insert(candidate.rows.end(), span.begin(), span.end());
            if (right_seen) append_part(right);
            make_recipes(candidate);
            for (auto recipe : candidate.recipes) {
                auto key = std::make_tuple(recipe.pattern, recipe.length, recipe.left_margin, recipe.right_margin);
                if (seen.insert(key).second) {
                    recipe.id = "candidate_barcode_" + std::to_string(++serial);
                    result.recipes.push_back(std::move(recipe));
                }
            }
        }
    }
    if (!result.recipes.empty()) result.notes.push_back("Barcode extraction uses known adapter-side recipes while the full layout remains unresolved");
}

} // namespace layout_detail

LayoutInference infer_layout(const std::vector<Read>& reads, const Options& options) {
    using namespace layout_detail;
    LayoutInference result;
    result.complete = true;
    auto vocab = vocabulary(options, result.notes);
    const int n = static_cast<int>(vocab.parts.size());
    std::vector<std::map<int, Hit>> observations;
    std::vector<std::size_t> support(n, 0), forward(n, 0), reverse(n, 0);
    std::vector<std::vector<double>> relative(n);
    std::vector<std::vector<std::size_t>> precedes(n, std::vector<std::size_t>(n, 0));
    std::map<std::pair<int, int>, std::vector<int>> gaps;
    std::vector<std::size_t> read_lengths;
    std::size_t repeated_reads = 0;
    for (const auto& read : reads) {
        if (options.boundary_poll) options.boundary_poll();
        auto hits = find_hits(read.sequence, vocab, options.error_ratio);
        bool repeated = false;
        for (std::size_t a = 0; a < hits.size() && !repeated; ++a) {
            if (vocab.parts[hits[a].part].poly) continue;
            for (std::size_t b = a + 1; b < hits.size(); ++b)
                if (hits[a].part == hits[b].part &&
                    (hits[a].end <= hits[b].start || hits[b].end <= hits[a].start)) {
                    repeated = true; break;
                }
        }
        repeated_reads += repeated;
        int f = 0, r = 0;
        for (const auto& hit : hits) if (!vocab.parts[hit.part].poly)
            (hit.reverse ? r : f) += static_cast<int>(vocab.parts[hit.part].seq.size()) - 3 * hit.edits;
        bool rc = r > f;
        std::map<int, Hit> selected;
        for (auto hit : hits) {
            if (hit.reverse != rc) continue;
            if (rc) { int old_start = hit.start; hit.start = read.sequence.size() - hit.end; hit.end = read.sequence.size() - old_start; }
            auto existing = selected.find(hit.part);
            if (existing == selected.end() || hit.start < existing->second.start) selected[hit.part] = hit;
        }
        for (const auto& item : selected) {
            int p = item.first;
            ++support[p]; ++(rc ? reverse[p] : forward[p]);
            relative[p].push_back(static_cast<double>(item.second.start) / std::max<std::size_t>(1, read.sequence.size()));
            for (const auto& other : selected) {
                if (item.second.end <= other.second.start) {
                    ++precedes[p][other.first];
                    gaps[{p, other.first}].push_back(other.second.start - item.second.end);
                }
            }
        }
        observations.push_back(std::move(selected)); read_lengths.push_back(read.sequence.size());
    }
    const std::size_t floor = std::max<std::size_t>(3, static_cast<std::size_t>(std::ceil(reads.size() * 0.10)));
    if (repeated_reads >= floor) {
        result.complete = false;
        result.notes.push_back("Repeated literal parts suggest concatemers or artifact units; a single backbone is unresolved");
    }
    std::vector<int> selected;
    for (int p = 0; p < n; ++p) {
        if (support[p]) result.parts.push_back({vocab.parts[p].id, vocab.parts[p].seq, support[p], forward[p], reverse[p]});
        if (support[p] >= floor) selected.push_back(p);
    }
    if (selected.empty()) {
        result.complete = false; result.notes.push_back("No known static parts reached the minimum read support");
        return result;
    }
    std::size_t literal_parts = 0;
    for (int part : selected) literal_parts += vocab.parts[part].poly == 0;
    if (literal_parts < 2) {
        result.complete = false;
        result.notes.push_back("At least two literal anchors are required for a resolved scaffold; poly runs alone cannot identify it");
    }
    // Pairwise same-read order supplies the backbone. Incompatible orders are explicit partial results.
    std::map<int, std::set<int>> edges;
    std::map<int, int> indegree;
    for (int p : selected) indegree[p] = 0;
    for (int a : selected) for (int b : selected) {
        if (a >= b) continue;
        auto ab = precedes[a][b], ba = precedes[b][a];
        if (ab + ba < 3) continue;
        int from = ab > ba ? a : b, to = ab > ba ? b : a;
        if (std::max(ab, ba) * 5 < (ab + ba) * 4) {
            result.complete = false;
            result.notes.push_back("Conflicting part orders suggest a mixture, repeats, or insufficient structure evidence");
            continue;
        }
        if (edges[from].insert(to).second) ++indegree[to];
    }
    std::vector<int> order;
    while (order.size() < selected.size()) {
        std::vector<int> ready;
        for (int p : selected) if (indegree[p] == 0 && std::find(order.begin(), order.end(), p) == order.end()) ready.push_back(p);
        if (ready.empty()) { result.complete = false; result.notes.push_back("Part-order graph contains an unresolved cycle"); break; }
        if (ready.size() > 1) { result.complete = false; result.notes.push_back("Some selected parts have no resolved relative order"); }
        std::sort(ready.begin(), ready.end(), [&](int a, int b) {
            auto avg = [&](int p) { return std::accumulate(relative[p].begin(), relative[p].end(), 0.0) / relative[p].size(); };
            return std::make_pair(avg(a), a) < std::make_pair(avg(b), b);
        });
        int next = ready.front(); order.push_back(next);
        for (int to : edges[next]) --indegree[to];
    }
    if (order.empty()) return result;
    if (order.size() < 2) { result.complete = false; result.notes.push_back("A single observed part does not establish a complete layout"); }
    std::vector<int> prefix, suffix;
    for (std::size_t i = 0; i < observations.size(); ++i) {
        auto first = observations[i].find(order.front()), last = observations[i].find(order.back());
        if (first != observations[i].end()) prefix.push_back(first->second.start);
        if (last != observations[i].end()) suffix.push_back(read_lengths[i] - last->second.end);
    }
    int serial = 0;
    auto append_gap = [&](int a, int b, const std::vector<int>& values) {
        auto rows = infer_gap(vocab, a, b, describe_gap(values), result, serial);
        result.rows.insert(result.rows.end(), rows.begin(), rows.end());
    };
    append_gap(-1, order.front(), prefix);
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto& part = vocab.parts[order[i]];
        Row row;
        row.id = part.id; row.sequence = part.poly ? "" : part.seq;
        row.type = "static"; row.klass = part.poly ? (part.poly == 'T' ? "poly_t" : "poly_a") : part.klass;
        row.length = part.poly ? 12 : part.seq.size(); row.length_spec = std::to_string(row.length);
        result.rows.push_back(std::move(row));
        if (i + 1 < order.size()) append_gap(order[i], order[i + 1], gaps[{order[i], order[i + 1]}]);
    }
    append_gap(order.back(), -2, suffix);
    std::set<std::string> ids;
    for (auto& row : result.rows) {
        std::string base = row.id.empty() ? "region" : row.id;
        row.id = base;
        for (int k = 2; !ids.insert(row.id).second; ++k) row.id = base + "_" + std::to_string(k);
    }
    make_recipes(result);
    fallback_recipes(vocab, support, result);
    bool has_insert = false;
    for (const auto& row : result.rows) has_insert |= row.klass == "read";
    if (!has_insert) {
        result.complete = false;
        result.notes.push_back("An extractable insert region has not been established");
    }
    result.notes.push_back("Barcode/UMI roles use known flanking-part contracts; catalog evidence is checked when --wl is enabled");
    std::sort(result.notes.begin(), result.notes.end());
    result.notes.erase(std::unique(result.notes.begin(), result.notes.end()), result.notes.end());
    result.notes.push_back("Known static parts supply the vocabulary; order and gap lengths are measured from this sample");
    return result;
}

void write_layout(const std::string& path, const LayoutInference& inference,
                  const std::vector<CatalogMatch>& matches) {
    using namespace layout_detail;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write discovered layout: " + path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    bool barcode = false;
    for (const auto& row : inference.rows) barcode |= row.klass == "barcode";
    out << (barcode ? "Read Layout\n" : "Read Layout:bulk\n")
        << "id,seq,expected_length,type,class,whitelist,flags\n";
    for (auto row : inference.rows) {
        const CatalogMatch* best = nullptr; bool tie = false;
        for (const auto& match : matches) {
            if (match.recipe_id != row.id || match.reads == 0) continue;
            if (!best || match.reads > best->reads) { best = &match; tie = false; }
            else if (match.reads == best->reads && match.path != best->path) tie = true;
        }
        if (best && !tie && !best->aliases.empty()) {
            row.whitelist = best->path;
            row.length = best->length; row.length_spec = std::to_string(best->length);
        }
        out << csv_field(row.id) << ',' << csv_field(row.sequence) << ','
            << csv_field(row.length_spec.empty() && row.length ? std::to_string(row.length) : row.length_spec) << ','
            << csv_field(row.type) << ',' << csv_field(row.klass) << ',' << csv_field(row.whitelist) << ',' << csv_field(row.flags) << '\n';
    }
    out.close();
}

} // namespace rad_discovery

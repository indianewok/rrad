#pragma once


struct spatial_reformat_options {
    std::string input, output, sample, map_path, bc1_path, bc2_path, metadata_path;
    std::string unit = "bin";
    int bin_size = 2;
    bool tissue_only = false;
    bool strict = false;
    std::function<void()> boundary_poll;
};

struct spatial_reformat_result {
    size_t records_read = 0, records_written = 0;
    std::map<std::string, size_t> counts;
};
spatial_reformat_result reformat_spatial(const spatial_reformat_options& options);

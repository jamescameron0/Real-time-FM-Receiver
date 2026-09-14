#include "RDSApplicationLayer.h"
#include "rds_helpers.h"

#include <iostream>
#include <vector>
#include <array>
#include <string>

// Helper function converting bitstream input to integer.
int bits_to_int(const std::vector<int>& bits) {
    int result = 0;
    for (int b : bits) {
        // result = (result << 1) | int(b)
        result = (result << 1) | b;
    }
    return result;
}

// Table defining PTY names
const std::array<std::string, 32> PTY_TABLE = {
    "None", "News", "Information", "Sports", "Talk", "Rock",
    "Classic Rock", "Adult Hits", "Soft Rock", "Top 40", "Country", "Oldies",
    "Soft music", "Nostalgia", "Jazz", "Classical", "Rhythm & Blues", "Soft R&B",
    "Language", "Religious Music", "Religious Talk", "Personality", "Public", "College",
    "Spanish Talk", "Spanish Music", "Hip Hop", "", "", "Weather", "Emergency Test", "Emergency"
};

// Default Constructor
RDSApplicationLayer::RDSApplicationLayer() {
    for (auto& c : ps_chars) c = '_';
}

// Converts input to hex code for Pi Code output
std::string RDSApplicationLayer::to_hex_4(int val) const {
    const char* digits = "0123456789ABCDEF";
    std::string out = "0000";
    out[0] = digits[(val >> 12) & 0xF];
    out[1] = digits[(val >> 8) & 0xF];
    out[2] = digits[(val >> 4) & 0xF];
    out[3] = digits[val & 0xF];
    return out;
}

// Extracts the 16 information bits from latter part of window. Returns window if less than 16 bits
std::vector<int> RDSApplicationLayer::info_bits(const std::vector<int>& window) {
    if (window.size() < 16) return window;
    return std::vector<int>(window.begin(), window.begin() + 16);
}

// Helper function for python slice function. Acts the same.
std::vector<int> RDSApplicationLayer::slice(const std::vector<int>& bits, int start, int end) {
    if (start >= bits.size()) return {};
    auto first = bits.begin() + start;
    auto last = bits.begin() + (end > (int)bits.size() ? bits.size() : end);
    return std::vector<int>(first, last);
}

// Takes in group block data and determines PI code, PTY value, and PS
void RDSApplicationLayer::decode_group(const std::vector<BlockData>& group) {
    // Extracts block bits from window
    std::vector<int> blk_a_bits = info_bits(group[0].window_bits);
    std::vector<int> blk_b_bits = info_bits(group[1].window_bits);
    std::vector<int> blk_d_bits = info_bits(group[3].window_bits);

    // Block A: PI
    pi_code = bits_to_int(blk_a_bits);

    // Block B - All bit ranges obtained in RBDS manual
    int group_type = bits_to_int(slice(blk_b_bits, 0, 4));
    int version    = blk_b_bits[4];
    int pty_val    = bits_to_int(slice(blk_b_bits, 6, 11));

    pty = pty_val;
    pty_name = PTY_TABLE[pty_val];

    // Only writes PS name if in Group 0 Version 0
    if (group_type == 0 && version == 0) {
        // information on PS name segment for current info bits stored in last bits of block B
        int seg = bits_to_int(slice(blk_b_bits, 14, 16));

        // PS info bits stored in block D
        int c1 = bits_to_int(slice(blk_d_bits, 0, 8));
        int c2 = bits_to_int(slice(blk_d_bits, 8, 16));

        // Stores first and second character from info bits into corresponding location based on segment
        ps_chars[seg * 2] = static_cast<char>(c1);
        ps_chars[seg * 2 + 1] = static_cast<char>(c2);

        // Writes all of ps_chars to string to be output
        ps_name = std::string(ps_chars.data(), 8);
    } 
}

// Adds block and block info to queue to be used in decoding
void RDSApplicationLayer::add_block(const std::vector<int>& window, RDSBlockOffset blk_offset) {
    // If first block, clear all prior block info and push block A
    if (blk_offset == A) {
        group_buf.clear();
        group_buf.push_back({window, blk_offset});
        return;
    }

    // Don't push info if queue empty
    if (group_buf.empty()) return;

    // If block B or C, push onto queue
    group_buf.push_back({window, blk_offset});

    // If block D, puch onto queue, then decode whole group is group is size 4. Then clear queue
    if (blk_offset == D && group_buf.size() == 4) {
        decode_group(group_buf);
        group_buf.clear();
    }
}

// Function to display RDS Info
void RDSApplicationLayer::display() {
    // New line to look nice
    std::cerr << "\033[2J\033[1;1H";

    // Converts PI code to hex
    std::string pi_str = (pi_code != 0) ? to_hex_4(pi_code) : "N/A";

    std::cerr << "  PI code: " << pi_str << "\n";
    std::cerr << "  Program type: " << pty << " – " << pty_name << "\n";
    std::cerr << "  Program service: " << ps_name << std::endl;
}
#ifndef DY4_RDSAPPLICATION_H
#define DY4_RDSAPPLICATION_H

#include <vector>
#include <string>
#include <array>

// Forward declaration or inclusion of rds_helpers if needed
// #include "rds_helpers.h" 

enum RDSBlockOffset { A, B, C, D };
struct BlockData {
    std::vector<int> window_bits; 
    RDSBlockOffset offset;
};

class RDSApplicationLayer {
private:
    int pi_code = 0; 
    int pty = 0;
    std::string pty_name = "N/A";
    
    std::array<char, 8> ps_chars;
    std::string ps_name = "________";
    
    std::vector<BlockData> group_buf;

    std::string to_hex_4(int val) const;
    std::vector<int> info_bits(const std::vector<int>& window);
    std::vector<int> slice(const std::vector<int>& bits, int start,  int end);
    void decode_group(const std::vector<BlockData>& group);

public:

    RDSApplicationLayer();

    void add_block(const std::vector<int>& window, RDSBlockOffset blk_offset);

    void display();
};

#endif
#ifndef RDS_FRAME_SYNCHRO_H
#define RDS_FRAME_SYNCHRO_H

#include <vector>
#include <array>
#include <map>

enum class Offset { SYNC, A, B, C, Cp, D };

class RDSFrameSynchro {
private:
    Offset state;
    
    using Syndrome = std::array<int, 10>;

    const std::map<Offset, Syndrome> SYNDROMES = {
        {Offset::A,  {1,1,1,1,0,1,1,0,0,0}},
        {Offset::B,  {1,1,1,1,0,1,0,1,0,0}},
        {Offset::C,  {1,0,0,1,0,1,1,1,0,0}},
        {Offset::Cp, {1,1,1,1,0,0,1,1,0,0}},
        {Offset::D,  {1,0,0,1,0,1,1,0,0,0}}
    };

    Syndrome compute_syndrome(const std::vector<int>& window);
    Offset match_syndrome(const Syndrome& syndrome);

public:
    RDSFrameSynchro();
    ~RDSFrameSynchro();

    std::pair<std::vector<std::pair<std::vector<int>, Offset>>, std::vector<int>> 
    process_block(std::vector<int>& bitstream);
};

#endif

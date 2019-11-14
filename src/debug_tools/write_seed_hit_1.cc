#include <debug_tools/write_seed_hit_1.h>

#include <fstream>

namespace raptor {

void WriteSeedHits(const std::string& out_path,
                   const std::vector<mindex::SeedHitPacked>& seed_hits, int32_t seed_len,
                   const std::string& qname, int64_t qlen,
                   const std::string& rname, int64_t rlen) {

    std::ofstream ofs(out_path);

    if (ofs.is_open() == false) {
        return;
    }

    ofs << qname.c_str() << "\t0\t" << qlen << "\t" << rname.c_str() << "\t0\t" << rlen << "\t0.0" << std::endl;

//   fprintf (fp, "%s\t0\t%ld\t%s\t0\t%ld\t0.0\n", qname.c_str(), qlen, rname.c_str(), rlen);

    for (size_t j=0; j<seed_hits.size(); j++) {
        int32_t cluster_id = seed_hits[j].TargetId() * 2 + (seed_hits[j].TargetRev() ? 1 : 0);
        ofs << seed_hits[j].QueryPos() << "\t" << seed_hits[j].TargetPos() << "\t" << cluster_id << std::endl;
        ofs << seed_hits[j].QueryPos() + seed_len << "\t" << seed_hits[j].TargetPos() + seed_len << "\t" << cluster_id << std::endl;
    }
}

void WriteTargetHits(const std::string& out_path,
                     const std::vector<raptor::ChainPtr>& target_hits,
                     int32_t seed_len,
                     const std::string& qname, int64_t qlen,
                     const std::string& rname, int64_t rlen) {

    std::ofstream ofs(out_path);

    if (ofs.is_open() == false) {
        return;
    }

    ofs << qname.c_str() << "\t0\t" << qlen << "\t" << rname.c_str() << "\t0\t" << rlen << "\t0.0" << std::endl;

    for (size_t i = 0; i < target_hits.size(); i++) {
        for (size_t j = 0; j < target_hits[i]->hits().size(); j++) {
            auto& hit = target_hits[i]->hits()[j];
            int32_t cluster_id = i;
            ofs << hit.QueryPos() << "\t" << hit.TargetPos() << "\t" << cluster_id << std::endl;
            ofs << hit.QueryPos() + seed_len << "\t" << hit.TargetPos() + seed_len << "\t" << cluster_id << std::endl;
        }
    }
}

}

#ifndef _IPA_INCLUDE_OVERLAP_FILE_H_
#define _IPA_INCLUDE_OVERLAP_FILE_H_

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <raptor_fetch/overlaps/overlap.h>
#include <raptor_fetch/overlaps/overlap_compact.h>
#include <raptor_fetch/overlaps/overlap_file_handlers.hpp>
#include <raptor_fetch/overlaps/overlap_parser.h>

namespace raptor {

class OverlapFile;
class QueryData;

using OverlapFilePtr = std::unique_ptr<raptor::OverlapFile>;
using QueryDataPtr = std::unique_ptr<raptor::QueryData>;

constexpr uint32_t FLAG_QUERY_LOW_LEN     = ((static_cast<uint32_t>(1)) << 0);  // 1
constexpr uint32_t FLAG_QUERY_CHIMERIC    = ((static_cast<uint32_t>(1)) << 1);  // 2
constexpr uint32_t FLAG_QUERY_CONTAINED   = ((static_cast<uint32_t>(1)) << 2);  // 4
constexpr uint32_t FLAG_QUERY_FILTERED    = ((static_cast<uint32_t>(1)) << 3);  // 8

OverlapFilePtr createOverlapFile(const std::string& path);

int64_t FindQuerySpan(const std::vector<raptor::OverlapCompactPtr>& overlaps, int64_t start);

class QueryData {
  public:
    QueryData()
        :  name(), len(0), clip_start(0), clip_end(0), flag(0), info()
    { }
    QueryData(const std::string& _name, int32_t _len, int32_t _clip_start, int32_t _clip_end, int32_t _flag)
        :  name(_name), len(_len), clip_start(_clip_start), clip_end(_clip_end), flag(_flag), info()
    { }

    std::string Verbose() const {
        std::ostringstream oss;
        oss << name << " " << len << " " << clip_start << " " << clip_end << " " << flag;
        return oss.str();
    }

    bool IsFiltered() const {
        return (flag != 0);
    }
    bool IsLowLen() const {
        return (flag & FLAG_QUERY_LOW_LEN);
    }
    bool IsChimeric() const {
        return (flag & FLAG_QUERY_CHIMERIC);
    }
    bool IsContained() const {
        return (flag & FLAG_QUERY_CONTAINED);
    }

    void SetFlagLowLen(bool val) {
        SetFlag_(FLAG_QUERY_LOW_LEN, val);
    }
    void SetFlagContained(bool val) {
        SetFlag_(FLAG_QUERY_CONTAINED, val);
    }
    void SetFlagFiltered(bool val) {
        SetFlag_(FLAG_QUERY_FILTERED, val);
    }
    void SetFlagChimeric(bool val) {
        SetFlag_(FLAG_QUERY_CHIMERIC, val);
    }

    std::string name;
    int32_t len;
    int32_t clip_start;
    int32_t clip_end;
    int32_t flag;
    std::string info;   // Additional info data can be stored here.

  private:
    void SetFlag_(uint32_t _flag, bool val) {
        if (val) {  flag |= _flag; }
        else {      flag &= ~_flag; }
    }
};

class OverlapFile {
  public:
    friend OverlapFilePtr createOverlapFile(const std::string& path);
    ~OverlapFile() = default;

    void OpenFile(const std::string& path);
    bool HashAllQnames();
    bool LoadNextBatch(size_t batch_size, bool single_query_only, int32_t min_len, double min_score, double min_idt);

    /*
     * Const getters.
    */
    const std::vector<std::unique_ptr<raptor::OverlapCompact>>& batch() const { return batch_; }
    const std::unordered_map<std::string, int32_t>& qname_map() const { return qname_map_; }
    const std::vector<std::unique_ptr<QueryData>>& query_data() const { return query_data_; }
    bool empty() { return batch_.empty(); }

    /*
     * Non-const getters.
    */
    std::vector<std::unique_ptr<QueryData>>& query_data() { return query_data_; }
    std::vector<std::unique_ptr<raptor::OverlapCompact>>& batch() { return batch_; }

    std::unique_ptr<QueryData>& GetQueryDataById(int32_t query_id) {
        size_t query_id_u = static_cast<size_t>(query_id);
        // std::cerr << "query_id = " << query_id << ", batch_start_id_ = " << batch_start_id_ << ", query_id_u = " << query_id_u << ", query_data_.size() = " << query_data_.size() << "\n";
        if (query_id < 0 || query_id_u >= query_data_.size()) {
            return dummy_null_query_data_;
        }
        return query_data_[query_id_u];
    }


  private:
    OverlapFile()
        :
            fp_handlers_(nullptr),
            qname_map_(),
            query_data_(),
            dummy_null_query_data_(nullptr),
            batch_(),
            batch_start_id_(0)
    { }

    std::unique_ptr<raptor::Overlap> YieldOverlap_(const std::string& line) const;

    raptor::OverlapFileHandlersPtr fp_handlers_;
    std::unordered_map<std::string, int32_t> qname_map_;
    std::vector<std::unique_ptr<QueryData>> query_data_;
    std::unique_ptr<QueryData> dummy_null_query_data_;

    std::vector<std::unique_ptr<raptor::OverlapCompact>> batch_;
    int32_t batch_start_id_;
};

}

#endif

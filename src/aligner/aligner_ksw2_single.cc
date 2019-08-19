#include "aligner/aligner_ksw2_single.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
// #include "ksw2/kseq.h"
#include "aligner/aligner_util.hpp"

#include <iostream>

// KSEQ_INIT(gzFile, gzread)

namespace raptor {

extern uint8_t seq_nt4_table[256];

//  = {
// 	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
// };

std::shared_ptr<AlignerBase> createAlignerKSW2Single(const raptor::AlignmentOptions &opt) {
    return std::shared_ptr<AlignerBase>(new AlignerKSW2Single(opt));
}

// static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
// {
// 	printf("%s\t%s\t%d", tname, qname, ez->score);
// 	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
// 	if (ez->n_cigar > 0) {
// 		int i;
// 		putchar('\t');
// 		for (i = 0; i < ez->n_cigar; ++i)
// 			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
// 	}
// 	putchar('\n');
// }

AlignerKSW2Single::AlignerKSW2Single(const raptor::AlignmentOptions &opt) : opt_(opt) {}

AlignerKSW2Single::~AlignerKSW2Single() {}

std::shared_ptr<raptor::AlignmentResult> AlignerKSW2Single::Global(const char *qseq, int64_t qlen,
                                                                   const char *tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::AlignerFailure);

    void *km = 0;
    ksw_extz_t ez;  // Alignment result.
    int flag = 0, zdrop = -1;

#ifdef HAVE_KALLOC
    km = km_init();
#endif

    memset(&ez, 0, sizeof(ksw_extz_t));

    auto mat = GenerateSimpleMatchMatrix<int8_t>((int8_t)opt_.p.match, (int8_t)opt_.p.mismatch, 5);
    // In Raptor definition, penalties are negative. KSW2 expects positive values.
    int8_t q = -opt_.p.w[0].open;  // Gap open. The intercept component of the affine function.
    int8_t e = -opt_.p.w[0].ext;  // Gap extend. The slope of the affine function.

    KSW2GlobalAlnWrapper_(km, (const int8_t *)qseq, qlen, (const int8_t *)tseq, tlen, 5, &mat[0], q,
                          e, opt_.bandwidth, zdrop, opt_.end_bonus, flag, &ez);

    result->score(ez.score);
    result->position(raptor::AlignmentPosition(0, qlen, 0, tlen));
    result->final_band(-1);
    result->status(raptor::AlignmentReturnValue::OK);

    std::vector<raptor::CigarOp> basic_cigar;
    for (size_t i = 0; i < ez.n_cigar; i++) {
        basic_cigar.push_back(raptor::CigarOp("MID"[ez.cigar[i] & 0xf], ez.cigar[i] >> 4));
    }
    result->cigar(raptor::ConvertBasicToExtCIGAR(qseq, qlen, tseq, tlen, basic_cigar));
    result->edit_dist(EditDistFromExtCIGAR(result->cigar()));

    if (result->score() == KSW_NEG_INF) {
        result->status(raptor::AlignmentReturnValue::Suboptimal);
    }

    kfree(km, ez.cigar);
#ifdef HAVE_KALLOC
    km_destroy(km);
#endif

    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerKSW2Single::Extend(const char *qseq, int64_t qlen,
                                                                   const char *tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    void *km = 0;
    ksw_extz_t ez;                                    // Alignment result.
    int flag = KSW_EZ_SCORE_ONLY | KSW_EZ_EXTZ_ONLY;  // | KSW_EZ_APPROX_DROP;
                                                      // int flag = KSW_EZ_SCORE_ONLY;

#ifdef HAVE_KALLOC
    km = km_init();
#endif

    memset(&ez, 0, sizeof(ksw_extz_t));

    auto mat = GenerateSimpleMatchMatrix<int8_t>((int8_t)opt_.p.match, (int8_t)opt_.p.mismatch, 5);
    // In Raptor definition, penalties are negative. KSW2 expects positive values for affine pieces.
    int8_t q = -opt_.p.w[0].open;  // Gap open. The intercept component of the affine function.
    int8_t e = -opt_.p.w[0].ext;  // Gap extend. The slope of the affine function.

    KSW2GlobalAlnWrapper_(km, (const int8_t *)qseq, qlen, (const int8_t *)tseq, tlen, 5, &mat[0], q,
                          e, opt_.zbandwidth, opt_.zdrop, opt_.end_bonus, flag, &ez);

    // print_aln("Query", "Target", &ez);

    result->score(ez.score);
    result->position(raptor::AlignmentPosition(0, ez.max_q + 1, 0, ez.max_t + 1));
    result->final_band(-1);
    result->status(raptor::AlignmentReturnValue::OK);
    result->max_score(ez.max);
    result->max_q_pos(ez.max_q);
    result->max_t_pos(ez.max_t);

    // std::cerr << "(Extend) mqe = " << ez.mqe << ", mqe_t = " << ez.mqe_t << ", mte = " << ez.mte
    // << ", mte_q = " << ez.mte_q
    // << ", ez.max_q = " << ez.max_q << ", ez.max_t = " << ez.max_t << std::endl;

    std::vector<raptor::CigarOp> basic_cigar;
    for (size_t i = 0; i < ez.n_cigar; i++) {
        basic_cigar.push_back(raptor::CigarOp("MID"[ez.cigar[i] & 0xf], ez.cigar[i] >> 4));
    }

    result->cigar(raptor::ConvertBasicToExtCIGAR(qseq, qlen, tseq, tlen, basic_cigar));
    result->edit_dist(EditDistFromExtCIGAR(result->cigar()));

    if (result->score() == KSW_NEG_INF) {
        result->status(raptor::AlignmentReturnValue::Suboptimal);
    }

    // printf ("Converted CIGAR:\n");
    // for (size_t i=0; i<result_->cigar.size(); i++) {
    //   printf ("%d%c", result_->cigar[i].count, result_->cigar[i].op);
    // }
    // printf ("\n");
    // printf ("Edit distance: %ld\n", result_->edit_dist);

    kfree(km, ez.cigar);
#ifdef HAVE_KALLOC
    km_destroy(km);
#endif

    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerKSW2Single::Local(const char *qseq, int64_t qlen,
                                                                  const char *tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::NotImplementedYet);
    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerKSW2Single::Semiglobal(const char *qseq,
                                                                       int64_t qlen,
                                                                       const char *tseq,
                                                                       int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::NotImplementedYet);
    return result;
}

void AlignerKSW2Single::KSW2GlobalAlnWrapper_(void *km, const int8_t *qseq_, int qlen,
                                              const int8_t *tseq_, int tlen, int8_t m,
                                              const int8_t *mat, int8_t q, int8_t e, int w,
                                              int zdrop, int end_bonus, int flag, ksw_extz_t *ez) {
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0;

    auto qseq = ConvertSeqAlphabet(qseq_, qlen, &seq_nt4_table[0]);
    auto tseq = ConvertSeqAlphabet(tseq_, tlen, &seq_nt4_table[0]);

    ksw_extz2_sse(km, qlen, (const uint8_t *)&qseq[0], tlen, (const uint8_t *)&tseq[0], m, mat, q,
                  e, w, zdrop, end_bonus, flag, ez);

    // const char *algo = "extd2_sse";
    // if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, (const uint8_t*)&qseq[0], tlen,
    // (const uint8_t*)&tseq[0], m, mat, q, e, w, zdrop, flag, ez); else if (strcmp(algo,
    // "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, (const uint8_t*)&qseq[0], tlen, (const
    // uint8_t*)&tseq[0], m, mat, q, e, q2, e2, w, zdrop, flag, ez);
    // // else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, (uint8_t*)qseq, tlen,
    // (uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez); else { 	fprintf(stderr, "ERROR: can't find
    // algorithm '%s'\n", algo); 	exit(1);
    // }
}

}  // namespace raptor

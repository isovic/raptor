/*
 * log_tools.h
 *
 *  Created on: Nov 10, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_LOG_LOG_TOOLS_H_
#define SRC_LOG_LOG_TOOLS_H_

#include <cstdarg>
#include <log/log_system.h>

#ifdef RAPTOR_TESTING_MODE
    #define DEBUG_QSEQ(params, qseq, log_call) \
        { \
            if (params->debug_qid == qseq->abs_id() || \
                params->debug_qname == std::string(qseq->header())) { \
                log_call; \
            } \
        }

    #define DEBUG_RUN(bool_do_debug, log_call) \
        { \
            if ((bool_do_debug)) { \
                log_call; \
            } \
        }

#else
    #define DEBUG_QSEQ(params, qseq, log_call) {}
    #define DEBUG_RUN(bool_do_debug, log_call) {}
#endif

#endif

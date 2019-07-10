/*
 * args_raptor.h
 *
 *  Created on: Mar 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_RAPTOR_H_
#define SRC_ARGS_RAPTOR_H_

#include <cstdint>
#include <memory>
#include <string>
#include <params/params_raptor.h>

namespace raptor {

void VerboseShortHelp(int argc, char **argv);
void VerboseShortHelpAndExit(int argc, char **argv, int ret_val = 1);
mindex::SequenceFormat WrapGetSequenceFormatFromPath(const mindex::SequenceFormat& apriori_in_fmt, const std::string& path);
std::vector<std::string> ExpandPathList(
                const mindex::SequenceFormat& apriori_in_fmt,
                const std::string& apriori_in_fmt_str,
                const std::vector<std::string>& in_paths);
void ValidateInputFiles(
                    int argc, char **argv,
                    const mindex::SequenceFormat& apriori_in_fmt,
                    const std::vector<std::string>& paths);
bool IsInputFormatRaptorDB(const mindex::SequenceFormat& apriori_in_fmt,
                            const std::vector<std::string>& paths);

int ProcessArgsRaptor(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptor> parameters);

} /* namespace raptor */

#endif /* SRC_ARGS_H_ */

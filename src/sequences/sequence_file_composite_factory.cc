/*
 * sequence_file_composite_factory.cc
 *
 *  Created on: Jul 22, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_composite_factory.h>
#include <sequences/sequence_file_enums.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>

namespace mindex {

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFactory(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt) {

	bool is_xml_used = false;

#ifdef RAPTOR_COMPILED_WITH_PBBAM
	// Sanity check for the input files. Allow only one XML file to be specified,
	// as by design of the SequenceFileCompositePbXml.
	for (const auto& in_path: in_paths) {
		if (mindex::GetSequenceFormatFromPath(in_path) == mindex::SequenceFormat::XML) {
			is_xml_used = true;
			break;
		}
	}
#endif

	mindex::SequenceFileCompositeBasePtr seq_file_parser = nullptr;

	// Set-up the parser for the correct sequence format.
	if (is_xml_used && in_paths.size() != 1) {
		FATAL_REPORT(ERR_UNEXPECTED_VALUE, "When using XML as input, only a single input can be specified.");

	} else if (is_xml_used && in_paths.size() == 1) {
        #ifdef RAPTOR_COMPILED_WITH_PBBAM
                seq_file_parser = mindex::createSequenceFileCompositePbXml(in_paths[0]);
        #else
                FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The raptor-reshape tool was not compiled with PacBio BAM support. Cannot use the .xml input files.");
        #endif

	} else {

		seq_file_parser = mindex::createSequenceFileCompositeFofn(in_paths, in_fmt);
	}

    return std::move(seq_file_parser);
}

}

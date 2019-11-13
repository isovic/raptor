/*
 * sequence_file_composite_factory.cc
 *
 *  Created on: Jul 22, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_composite_factory.h>
#include <sequences/sequence_file_enums.h>
#include <sequences/sequence_file_utils.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>

namespace mindex {

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFactory(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt) {
	// This checks that there is at most one RaptorDB or one XML file specified,
	// and that all input files actually exist on disk.
	bool rv = mindex::ValidateInputFiles(in_fmt, in_paths);

	if (rv == false) {
		FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Invalid input files.");
	}

	mindex::SequenceFileCompositeBasePtr seq_file_parser = nullptr;

	bool is_xml_used = false;
	#ifdef RAPTOR_COMPILED_WITH_PBBAM
		if (in_fmt == mindex::SequenceFormat::XML) {
				is_xml_used = true;
				seq_file_parser = mindex::createSequenceFileCompositePbXml(in_paths[0]);
		}
	#endif
	if (is_xml_used == false) {
		seq_file_parser = mindex::createSequenceFileCompositeFofn(in_paths, in_fmt);
	}

    return std::move(seq_file_parser);
}

}

/*
 * params_raptor_repack.cc
 *
 *  Created on: Dec 22, 2018
 *      Author: Ivan Sovic
 */

#include <raptor_reshape/params_raptor_reshape.h>

namespace raptor {

std::shared_ptr<raptor::ParamsRaptorReshape> createParamsRaptorReshape() {
    return std::shared_ptr<raptor::ParamsRaptorReshape>(new raptor::ParamsRaptorReshape);
}

ParamsRaptorReshape::ParamsRaptorReshape()
    : subprogram(),
      command_line(),
      verbose_level(0),
      in_paths(),
      out_prefix(),
      in_batch_size(400.0),
      keep_lowercase(false),
      rename_seqs(false),
      block_size(400.0),
      split_blocks(false),
      symlink_files(false)

{}

}  // namespace raptor

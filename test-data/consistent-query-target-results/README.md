This test data is intended to verify that both the A -> B and B -> A alignment produces the same identical result.
There was a bug which caused the different direction to have wrong alignment extension, due to the way how flank sequences were selected. The flanks were chosen with respect to query length.
Now, the flanks are chosen with respect to the shortes one (either query or target).


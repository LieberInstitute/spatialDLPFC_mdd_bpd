For *whatever* reason, `nnSVG` will complete with no errors for sample `V13B23-279_A1` if and only if(?) it is run with `verbose=T` and after another sample that does complete (`V13B23-279_C1`).<break><break>
I didn't actually test if both verbose and running after another sample are necessary but the proof that one of those two is required is laid out in the following logs.<break>
- "nnSVG_per-sample_RERUN_13764827_4294967294.log": Produces same BiocParallel error that occurred the first time `nnSVG` was run
- "nnSVG_BiocParallel-error_13764884.log": The only difference to the above run that errored is the use of `verbose=T` and running after another sample
<break><break>
Use this as a template to complete the `nnSVG` runs for the three remaining samples outlined in `processed-data/04_feature_selection/per-sample_spe_RERUN_list.txt`.

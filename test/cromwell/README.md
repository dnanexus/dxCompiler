# Excluded Tests

The following tests from the Cromwell suite will not pass on DNAnexus due to the following know limitations:

* Overriding task-level inputs is not supported
  * local_backend
  * string_interpolation
  * call_cache_hit_prefixes
  * declarations
  * reference_disk_test
  * optional_parameter
* Cannot reuse task/workflow names across imports
  * sub
  * sub_sub
  * echo
  * sub_workflow_no_output
  * recursive_imports
* No support for GCP (`gs://`) URIs
  * large_final_workflow_outputs_dir
  * input_from_bucket_with_requester_pays
  * input_expressions
* No support for the alternative workflow output syntax that has been deprecated since draft2
  * optional_declarations
  * sub_workflow_interactions
  * unscattered
* Forward references not supported
  * inter_scatter_dependencies

The following tests from the Cromwell suite are invalid and have been updated:

* Use `docker` keyword in runtime, but this is not allowed in version `development` (changed to `container`)
  * string_interpolation_optional
  * none_literal
* Assigning a String value to an Int variable
  * sub_workflow_one_output_import
  * sub_workflow_var_refs
  * sub_workflow_var_refs_import
* Reuse of names within the same scope
  * globbingBehavior

The following tests from the Cromwell suite are invalid and will not pass on DNAnexus for the following reasons:

* The write_tsv function does not support an `Array[Object]` argument
  * read_write_functions
* Version `development` does not allow arbitrary runtime keys:
  * afters_and_ifs
  * afters
  * afters_and_scatters
  * custom_cacheworthy_attributes

The following tests from the Cromwell suite require non-standard, Cromwell-specific behavior and will not pass on DNAnexus for the following reasons:

* Non-standard runtime key `backend` is not supported:
    * docker_alpine
    * parallel_composite_uploads_lib
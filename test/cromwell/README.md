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
* No support for GCP
  * large_final_workflow_outputs_dir
  * input_from_bucket_with_requester_pays
  * input_expressions
  * missing_delete
  * confirm_preemptible
  * call_cache_capoeira_jes
  * dedup_localizations_papi_v2
  * papi_v2_log
  * papi_v2_plain_detritus
  * monitoring_log
  * docker_size_dockerhub
  * docker_size_gcr
* No support for the alternative workflow output syntax that has been deprecated since draft2
  * optional_declarations
  * sub_workflow_interactions
  * unscattered
* Forward references not supported
  * inter_scatter_dependencies
* Cannot have workflows with no stages
  * recursive_imports_no_subwf

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
* DNAnexus does not support `continueOnReturnCode` - workflow updated to use standard `returnCodes` runtime attribute;
  * continue_on_return_code
  * exit

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
  * call_cache_capoeira_local
  * backendWithNoDocker
  * docker_image_cache_true
  * dummy_scatter
  * fofn_caching
  * hello_private_repo
  * local_bourne
  * papi_v2_gcsa
  * parallel_composite_uploads
  * call_cache_capoeira_tes
  * check_network_in_vpc
  * workbench_health_monitor_check
* Non-standard `volatile` meta key
  * monitoring_image_script
* TMPDIR is not a default environment variable in DNAnexus worker
  * tmp_dir
* DNAnexus puts a limit of 32k on string length
  * long_cmd
* Custom mount points are not supported
  * custom_mount_point
* Divide-by-zero is not short-circuited
  * short_circuit
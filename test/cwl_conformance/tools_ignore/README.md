These are tool tests that we are currently ignoring.

## Tests that can't be run due to bugs in cwltool or cwljava

* Enums in packed workflows are parsed to create URIs ([issue](https://github.com/common-workflow-language/cwltool/issues/1513))
    * anon_enum_inside_array
    * anon_enum_inside_array_inside_schemadef
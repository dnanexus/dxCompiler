# Adding a New DNAnexus Region

If dxCompiler needs to be enabled in a new DNAnexus region the following should be done in the staging & production environments:

* dxCompiler requires WDL and CWL assets to be stored in a public "dxCompiler\_region" project. For existing regions 
supported by DNAnexus - consult the internal Confluence page.
```bash
dx new project <PROJECT_NAME> --region <REGION> --bill-to=org-dnanexus_apps
```
* Make the project public
```bash
dx api project-xxxx invite '{"invitee": "PUBLIC", "level": "VIEW"}'
```
* The [release script](/scripts/build_release.py#L33) that tests and releases dxCompiler in all regions needs to be 
updated (one line change).
* The app used for copying the assets to different regions during a release ([app-dxwdl_copy](/scripts/dxcompiler_copy)) 
needs to be enabled in the new region (please update `regionalOptions`, `whatsNew`, and increment the `version` of the app).
  * To publish `dxwdl_copy` app (`Copy file`) go to `Actions` > `Build & publish 'dxCompiler Copy File' app` and click `Run workflow` on the right side.
* Update the script for [multi-region testing](/scripts/multi_region_tests.py#L24).

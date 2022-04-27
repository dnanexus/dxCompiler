

# Config file-based custom reorg applet example use case


Consider the following WDL task and workflow:

```
version 1.0

workflow simple_wf {
    input {}
    call simple_task
    output {
      Array[File] out = simple_task.out
    }
}

task simple_task {

  input {}

  command <<<
    echo "vcf out" > out.vcf
    echo "bam out" > out.bam
  >>>

  output {
     Array[File] out= ["out.vcf", "out.bam"]
  }
}
```

The workflow will generate an array of files - `out.vcf` and `out.bam`.
Suppose we want to reorganize the results into to 2 different folders on the project.

`out.vcf` will go to /folder_for_vcf
`out.bam` will go to /folder_for_bam

A simple reorg applet is shown below.

## `___reorg_conf` as conf.json file.

The key is the suffix for the files that we will move to the destination declared in the value of the JSON object.

```
{
  "bam": "/folder_for_bam",
  "vcf": "/folder_for_vcf"
}     
```

## code.py

The applet code does the following.

1) Read the configuration file.
2) Retrieve the list of output using the job and analysis id.
3) Move the files using dxpy.

```
import json
import dxpy
import time

@dxpy.entry_point('main')
def main(reorg_conf___=None, reorg_status___=None):

    # download and parse `reorg_conf___`
    conf_file = dxpy.DXFile(reorg_conf___)
    dxpy.download_dxfile(conf_file.get_id(), "conf.json")
    with open('conf.json') as f:
        conf = json.load(f)

    # find the output stage of the current analysis
    analysis_id = dxpy.describe(dxpy.JOB_ID)["analysis"]
    # describe the analysis in a loop until dependsOn is empty
    # or contains only this reorg job's ID
    while True:
        analysis_desc = dxpy.describe(analysis_id)
        depends_on = analysis_desc.get("dependsOn")
        if not depends_on or depends_on == [dxpy.JOB_ID]:
            break
        else:
            time.sleep(3)

    stages = analysis_desc["stages"]

    # retrieve the dictionary containing outputs, where key is the name of output and value is the link to the file.
    output_map = [x['execution']['output'] for x in stages if x['id'] == 'stage-outputs'][0]
    out = output_map['out']
    filenames = [dxpy.DXFile(x).describe(fields={'name'}) for x in out]
    bam = [x['id'] for x in filenames if x["name"].endswith('.bam')][0]
    vcf = [x['id'] for x in filenames if x["name"].endswith('.vcf')][0]

    vcf_folder = conf['vcf']
    bam_folder = conf['bam']

    # get the container instance
    dx_container = dxpy.DXProject(dxpy.PROJECT_CONTEXT_ID)
    dx_container.new_folder(vcf_folder, parents=True)
    dx_container.new_folder(bam_folder, parents=True)
    dx_container.move(
        destination=vcf_folder,
        objects=[ vcf ]
    )
    dx_container.move(
        destination=bam_folder,
        objects=[ bam ],
    )

    return dict(outputs=out)
```

## dxapp.json

This is the spec for the applet.
`reorg_conf___` is declared in `custom-reorg.conf` in extras.json (JSON file provided to `-extras`).

```
{
  "name": "custom_reorg_app",
  "title": "Example custom reorg app",
  "summary": "A small example to show how to use the config file based custom reorg app",
  "dxapi": "1.0.0",
  "inputSpec": [
    {
      "name": "reorg_conf___",
      "label": "Config",
      "help": "",
      "class": "file",
      "patterns": ["*"],
      "optional": true
    },
    {
      "name": "reorg_status___",
      "label": "Config",
      "help": "",
      "class": "string",
      "optional": true
    }
  ],
  "outputSpec": [
    {
      "name": "outputs",
      "label": "Outputs",
      "help": "",
      "class": "array:file",
      "patterns": ["*"]
    }
  ],
  "runSpec": {
    "version": "0",
    "interpreter": "python3",
    "timeoutPolicy": {
      "*": {
        "hours": 48
      }
    },
    "distribution": "Ubuntu",
    "release": "20.04",
    "file": "reorg.py"
  },
  "access": {
    "network": [
      "*"
    ],
    "project": "CONTRIBUTE"
  },
  "ignoreReuse": false,
  "regionalOptions": {
    "aws:us-east-1": {
      "systemRequirements": {
        "*": {
          "instanceType": "mem1_ssd1_v2_x4"
        }
      }
    }
  }
}
```

## `extras.json` (JSON file provided to `-extras`).

```
{
  "custom-reorg" : {
    "app_id" : "applet-12345678910",
    "conf" : "dx://file-xxxxxxxx"
  }
}

```

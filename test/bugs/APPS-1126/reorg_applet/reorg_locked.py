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
    print(stages)
    output_map = [x['execution']['output'] for x in stages if x['id'] == 'stage-0'][0]
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


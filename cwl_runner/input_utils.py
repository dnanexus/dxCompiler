import os
import sys
import json
import cache
import utils
import constants as C

def get_modified_input(i):
    if type(i) is dict: # File and Directory are always dicts. # TODO: is there a different class that would be dict as well?
        location = i.get("location", i.get("path"))
        cached_location = cache.get_cache(location)
        if cached_location is not None:
            i["location"] = cached_location
            i.pop("path", None)
            return i

        file_class = i.get("class")
        if location is not None and file_class is not None and not location.startswith("dx://"):
            if os.path.isdir(location) and file_class == "Directory":
                print("Directories are not supported!", file=sys.stderr)
                exit(33)
            elif os.path.isfile(location) and file_class == "File":
                platform_location = utils.upload_file(location, C.test_dx_folder)
                cache.add_to_cache(location, platform_location)
                i["location"] = platform_location
                i.pop("path", None)
            elif file_class in ["Directory", "File"]:
                raise FileNotFoundError(f"File/Folder {location} does not exists, are you in correct local folder?")
    return i

def create_dx_input(jobfile):
    # issues:
    #         if file does not exist, exception is thrown and no json is generated, even if some files were uploaded
    with open(jobfile) as input_file:
        inputs = json.load(input_file)
        for i in inputs:
            if type(inputs[i]) is list:
                for idx, value in enumerate(inputs[i]):
                    inputs[i][idx] = get_modified_input(value)
            modified_input = get_modified_input(inputs[i])
            inputs[i] = modified_input
    return inputs


def write_dx_input(new_input, jobfile):
    with open(jobfile, 'w') as dx_input:
        dx_input.write(json.dumps(new_input, indent=4))
    return jobfile


def get_new_dx_input(input_path):
    basename = os.path.basename(input_path)
    new_input_filename = f"{os.path.splitext(basename)[0]}.dx.json"
    return new_input_filename
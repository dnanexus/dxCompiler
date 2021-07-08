#!/usr/bin/python
# helper file for cromwell_utils.sh
import sys, re, json

regex = 'metadata\s*\{[\n\sa-zA-Z:_".?!@#$%^&*()]*}'
file_results = {}
test_file = sys.argv[1]
output_file = sys.argv[2]

with open(test_file) as file:
    data = file.read()
    m = re.search(regex, data)
    if m:
        metadata = m.group(0).splitlines()
        for line in metadata:
            line = line.strip()
            if line.startswith('"outputs'):
                regex_output = '"outputs\.(.*)":\s*(.*)'
                m = re.search(regex_output, line)
                if m:
                    result = m.group(2).strip()
                    if result.startswith('"') and result.endswith('"'):
                        result = result[1:-1]
                    elif result == "false":
                        result = False
                    elif result == "true":
                        result = True
                    elif result == "null":
                        result = None
                    elif "." in result:
                        try:
                            result = float(result)
                        except:
                            print(f"Error, please check the test file {test_file}. Using output string as an output type.")
                    else:
                        try:
                            result = int(result)
                        except:
                            print(f"Error, please check the test file {test_file}. Using output string as an output type.")
                    file_results[m.group(1)] = result

if file_results:
    json_obj = json.dumps(file_results, indent=4)
    folder, file_name = output_file.rsplit("/", 1)
    file_name = file_name[:-4] + "_results.json"
    with open(folder + "/" + file_name, "w") as results:
        results.write(json_obj)

""" Converts a MIBiG database dump in CSV form to a directory of entry JSONs """

import json
import os
import sys


def convert(csv_path: str, output_path: str) -> int:
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    with open(csv_path) as handle:
        lines = handle.readlines()
    for line in lines:
        data = json.loads(line.split("\t", 4)[-1])
        accession = data["general_params"]["mibig_accession"].strip().strip("'\"")
        with open(os.path.join(output_path, f"{accession}.json"), "w") as handle:
            json.dump(data, handle, indent=1)


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(f"Usage: {sys.argv[0]} output_dir [csv_file: defaults to 'mibig.csv']")
        sys.exit(1)
    try:
        convert(sys.argv[2] if len(sys.argv) > 2 else "mibig.csv", sys.argv[1])
    except PermissionError as err:
        print(err, file=sys.stderr)
        sys.exit(1)

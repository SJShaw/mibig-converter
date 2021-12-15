import glob
import json
import os
import sys

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} destination_mask source_dir")
    print(f"E.g.: {sys.argv[0]} 'outputs/*.json' existing/data/")
    sys.exit(1)

for dest in glob.glob(sys.argv[1]):
    base = os.path.basename(dest)
    with open(os.path.join(sys.argv[2], base)) as handle:
        source = json.load(handle)

    with open(dest) as handle:
        data = json.load(handle)

    data["cluster"]["ncbi_tax_id"] = source["cluster"]["ncbi_tax_id"]
    data["cluster"]["organism_name"] = source["cluster"]["organism_name"]
    for loci in data["cluster"]["loci"]:
        try:
            if loci["accession"].startswith(source["cluster"]["loci"]["accession"].split(".")[0]):
                loci["accession"] = source["cluster"]["loci"]["accession"]
        except TypeError:
            print(dest, source["cluster"]["loci"]["accession"].split(".")[0])
            print(dest, repr(loci))
            raise
    change_addition = data["changelog"]
    data["changelog"] = source["changelog"]
    if not any(change["version"] == "2.1" for change in source["changelog"]):
        data["changelog"].extend(change_addition)

    with open(dest, "w") as handle:
        json.dump(data, handle, indent=4, sort_keys=True, ensure_ascii=False)

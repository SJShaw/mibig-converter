import glob
import json
import os
from typing import Any, Dict, List, Union
import warnings

import jsonschema

JSON = Dict[str, Any]


class ConversionError(ValueError):
    pass


def _gather_schemas(schema_dir: str) -> Dict[str, None]:
    schemas = {}
    print(schema_dir)
    for schema_file in glob.glob(os.path.join(schema_dir, "*.json")):
        with open(schema_file) as handle:
            schema = json.load(handle)
        version = schema["$schema_version"]
        assert version not in schemas
        schemas[version] = schema
    return schemas


SCHEMAS = _gather_schemas(os.path.join(os.path.dirname(__file__), "schemas"))
LATEST = sorted(SCHEMAS)[-1]


def convert_single(input_path: str, output_path: str, target_schema: str = LATEST, mibig_version: str = "2.1") -> None:
    with open(input_path) as handle:
        before = json.load(handle)

    if "personal" in before:
        warnings.warn("dropping submitter information")  # TODO
        before.pop("personal")

    # drop things that aren't kept in any version
    for key in ["version"]:
        before.pop(key, None)

    # convert
    after = transform(before, mibig_version)

    # write
    with open(output_path, "w") as handle:
        json.dump(after, handle, indent=1)

    # ensure the result is valid according to the chosen schema
    try:
        is_valid(after, target_schema)
    except ValueError as err:
        accession = os.path.basename(input_path).rsplit(".", 1)[0]
        raise ConversionError(f"failure converting {accession}: {err}")


def convert(input_paths: Union[str, List[str]], output_dir: str) -> None:
    assert os.path.exists(output_dir) and os.path.isdir(output_dir)
    if isinstance(input_paths, str):
        return convert_single(input_paths, os.path.join(output_dir, os.path.basename(input_paths)))
    for i, path in enumerate(sorted(input_paths)):
        try:
            convert_single(path, os.path.join(output_dir, os.path.basename(path)))
        except (AssertionError, AttributeError, KeyError, ValueError):
            print(f"{i}: {os.path.basename(path)}")
            raise


def build_cluster(data: JSON) -> JSON:
    def convert_loci(cluster: JSON) -> JSON:
        old = cluster.pop("loci")
        assert len(old["nucl_acc"]) == 1
        accession_info = old["nucl_acc"][0]

        loci = {
            "completeness": old.pop("complete", "Unknown"),
            "accession": accession_info["Accession"],
            "start_coord": accession_info["start_coord"],
            "end_coord": accession_info["end_coord"],
        }
        warnings.warn("mixs_compliant not checked or added")  # TODO mixs_compliant -> bool

        evidence = accession_info.pop("conn_comp_cluster", [])
        if evidence:
            loci["evidence"] = evidence

        return loci

    cluster = data.pop("general_params")
    assert "mibig_accession" in cluster
    # mandatory members:
    # loci
    cluster["loci"] = convert_loci(cluster)
    # compounds
    cluster["compounds"] = convert_compounds(cluster.pop("compounds"))
    # publications
    publications = []
    existing = cluster["publications"]
    if isinstance(existing, str):
        existing = [ex.strip() for ex in existing.split(",")]
    for publication in existing:
        if publication.isdigit():
            publications.append(f"pubmed:{publication}")
        else:
            publications.append(publication)
    cluster["publications"] = publications
    # organism_name
    warnings.warn("missing organism name lookups")
    cluster["organism_name"] = "missing"  # TODO
    # ncbi_tax_id
    warnings.warn("missing NCBI tax ID lookups")
    cluster["ncbi_tax_id"] = "-1"  # TODO
    # minimal
    cluster["minimal"] = cluster.get("loci", {}).get("evidence") is None

    # optionals:
    if "NRP" in cluster:
        cluster["nrp"] = convert_NRP(cluster.pop("NRP"))
    # alkaloid
    # polyketide
    # other
    # loci
    if "genes" in cluster:
        cluster["genes"] = convert_genes(cluster.pop("genes"))
    # ripp
    # saccharide
    # terpene

    return cluster


def convert_compounds(old: JSON) -> JSON:
    def convert_compound(old_compound: JSON) -> JSON:
        new = {
            "compound": old_compound.pop("compound"),
        }
        optional_renames = [
            ("chem_synonyms", "chem_synonyms"),
            ("chem_acts", "chem_act"),
            ("chem_target", "chem_targets"),
            ("mol_mass", "mol_mass"),
            ("mass_ion_type", "mass_spec_ion_type"),
            ("chem_targets", "chem_target"),
            ("molecular_formula", "molecular_formula"),
            ("chem_struct", "chem_struct"),
            ("evidence", "evidence"),
            ("chem_moieties", "chem_moieties"),
            ("chem_act", "chem_acts"),
        ]
        for src, dest in optional_renames:
            if src in old_compound:
                new[dest] = old_compound.pop(src)

        dbs = []
        for db in old_compound.pop("databases_deposited", []):
            name = db.lower()
            old_id = old_compound.pop(f"{name}_id")
            dbs.append(f"{name}:{old_id}")
        if "chemspider_id" in old_compound:
            dbs.append(f"chemspider:{old_compound.pop('chemspider_id')}")
        if dbs:
            new["database_id"] = dbs
        old_compound.pop("database_deposited", None)

        assert not old_compound, old_compound
        return new

    compounds = []
    for compound in old:
        compounds.append(convert_compound(compound))
        assert not compound, compound
    return compounds


def convert_genes(old: JSON) -> JSON:
    def convert_gene(old_gene: JSON) -> JSON:
        return old_gene
    new = {}
    genes = []
    for gene in old.pop("gene", []):
        genes.append(convert_gene(gene))
    if genes:
        new["annotations"] = genes
    if "operon" in old:
        operons = old.pop("operon")
        new["operons"] = operons

    assert not old, old
    return new


def convert_NRP(old: JSON) -> JSON:
    new = {}

    if "lin_cycl_nrp" in old:
        new["cyclic"] = old.pop("lin_cycl_nrp") == "Cyclic"
    if "nrps_release_type" in old:
        new["release_type"] = [old.pop("nrps_release_type")]
    if "subclass" in old:
        new["subclass"] = old.pop("subclass")
    if "lipid_moiety" in old:
        new["lipid_moiety"] = old.pop("lipid_moiety")

    thios = []
    thio_type = old.pop("nrps_te_type")
    for thio in old.pop("nrps_thioesterase", []):
        thios.append({
            "gene": thio,
            "thioesterase_type": thio_type,
        })
    if thios:
        new["thioesterases"] = thios

    genes = []
    for gene in old.pop("nrps_genes", []):
        modules = gene.pop("modules", [])
        if modules:
            gene["modules"] = sorted(modules, key=lambda module: module["module_nr"])
        genes.append(gene)
    if genes:
        new["nrps_genes"] = genes

    assert not old, old
    return new


def transform(data: JSON, mibig_version: str) -> JSON:
    transformed = copy.deepcopy(data)

    transformed["cluster"] = build_cluster(transformed)

    changelog = [
        # TODO add submitter
        {
            "comments": [
                "Migrated from v1.4"
            ],
            "contributors": [
                "AAAAAAAAAAAAAAAAAAAAAAAA"
            ],
            "version": mibig_version
        }
    ]
    transformed["changelog"] = changelog


    return transformed


def is_valid(data: JSON, schema_version: str) -> bool:
    try:
        jsonschema.validate(data, SCHEMAS[schema_version])
    except jsonschema.ValidationError as err:
        summary = str(err).splitlines()[0]
        path = "][".join(map(str, list(err.path)))
        if path:
            raise ValueError(f"invalid result: {summary}: [{path}]")
        raise ValueError(f"invalid result: {summary}")
    return True

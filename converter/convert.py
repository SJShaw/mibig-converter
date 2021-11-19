import glob
import json
import os
from typing import Any, Dict, List, Tuple, Union
import warnings

import jsonschema

JSON = Dict[str, Any]

VALID_CYCLICS = {"Cyclic", "Unknown", "Linear"}


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
        json.dump(after, handle, indent=4, sort_keys=True)

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
            print(os.path.basename(path))
            convert_single(path, os.path.join(output_dir, os.path.basename(path)))
        except (AssertionError, AttributeError, KeyError, ValueError):
            print(f"{i}: {os.path.basename(path)}")
            raise


def pub_cmp_parts(pub):
    """ A comparison function for sorting publications by:
        - database, alphabetically; then
        - publication number, numerically and reversed

    """
    tag, value = pub.split(":", 1)
    if value.isdigit():
        value = int(value) * -1  # to reverse on the int portion only, so newer publications from the same source are at the top
    return tag, value


def build_cluster(old: JSON) -> JSON:
    def convert_loci(cluster: JSON) -> JSON:
        old = cluster.pop("loci")
        if len(old["nucl_acc"]) > 1:
            warnings.warn("multiple reference accessions detected, using first")
        accession_info = old.pop("nucl_acc")[0]

        completeness = old.pop("complete", "Unknown")
        comp_mapping = {
            "unknown": "Unknown",
            "partial": "incomplete",
        }
        completeness = comp_mapping.get(completeness, completeness)
        loci = {
            "completeness": completeness,
            "accession": accession_info["Accession"].strip(),
            "start_coord": accession_info["start_coord"],
            "end_coord": accession_info["end_coord"],
        }
        for coord in ["start_coord", "end_coord"]:
            if loci[coord] == -1:
                loci.pop(coord)
        if loci.get("start_coord") == 0:
            loci["start_coord"] = 1
        warnings.warn("mixs_compliant not checked or added")  # TODO mixs_compliant -> bool

        evidence = accession_info.pop("conn_comp_cluster", [])
        if evidence:
            if isinstance(evidence, str):
                evidence = [ev.strip() for ev in evidence.strip().split(",")]
            conversions = {
                "Proven expression in natural host": "Gene expression correlated with compound production",
                "Other": None,
            }
            new_evs = []
            for ev in evidence:
                new_ev = conversions.get(ev, ev)
                if new_ev is not None:
                    new_evs.append(new_ev)
            if new_evs:
                loci["evidence"] = new_evs

        assert not old, old
        return loci

    new = {}
    assert "mibig_accession" in old
    # mandatory members:
    # loci
    new["mibig_accession"] = old.pop("mibig_accession")
    new["biosyn_class"] = old.pop("biosyn_class")
    if "Nucleoside" in new["biosyn_class"]:
        new["biosyn_class"] = [cl for cl in new["biosyn_class"] if cl != "Nucleoside"]
        assert "Other" not in old or old["Other"]["other_subclass"] == "None"
        new["biosyn_class"].append("Other")
        old["Other"] = {"other_subclass": "nucleoside"}

    new["loci"] = convert_loci(old)
    # compounds
    new["compounds"] = convert_compounds(old.pop("compounds"))
    # publications
    publications = []
    existing = old.pop("publications", [])
    # some are strings, comma separated, some are colon separated, some are just strings
    if isinstance(existing, str):
        if "," in existing:
            existing = [ex.strip() for ex in existing.strip(" ,").split(",")]
        elif ";" in existing:
            existing = [ex.strip() for ex in existing.strip(" ;").split(";")]
        else:
            existing = [existing.strip()]
    # some are even lists of comma separated strings
    if isinstance(existing, list):
        new_pubs = []
        for pub in existing:
            if "," in pub:
                new_pubs.extend(ex for ex in pub.strip(" ,").split(","))
            elif ";" in pub:
                new_pubs.extend(ex.strip() for ex in pub.strip(" ;").split(";"))
            else:
                new_pubs.append(pub)
        existing = [pub.strip("' ") for pub in new_pubs]
    # some are buried in genes
    if "genes" in old:
        new_genes, gene_publications = convert_genes(old.pop("genes"))
        existing.extend(gene_publications)
        if new_genes:
            new["genes"] = new_genes
    # but eventually they should just be single strings in a list
    for publication in existing:
        assert "," not in publication, publication
        if publication.isdigit():
            publications.append(f"pubmed:{publication}")
        elif publication.startswith("PMC") and publication[3:].isdigit():
            publications.append(f"pubmed:{publication[3:]}")
        elif publication.startswith("10."):
            publications.append(f"doi:{publication}")
        elif publication.startswith("doi:"):
            publications.append(f"doi:{publication.split(':', 1)[1].strip()}")
        elif publication.startswith("doi.org/"):
            publications.append(f"doi:{publication.split('/', 1)[1].strip()}")
        elif publication.startswith("https://doi.org/"):
            publications.append(f"doi:{publication.split('/', 3)[3].strip()}")
        elif publication in ["-", "unpublished"]:
            continue
        else:
            publications.append(publication)
    if not publications:
        warnings.warn(f"{new['mibig_accession']} is missing publications, setting dummy")  # TODO
        publications.append("pubmed:000000")
    new["publications"] = sorted(publications, key=pub_cmp_parts)
    # organism_name
    warnings.warn("missing organism name lookups")
    new["organism_name"] = "missing"  # TODO
    # ncbi_tax_id
    warnings.warn("missing NCBI tax ID lookups")
    new["ncbi_tax_id"] = "-1"  # TODO
    # minimal
    new["minimal"] = old.pop("minimal", new.get("loci", {}).get("evidence") is None)

    # optionals:
    if "NRP" in old:
        new["nrp"] = convert_NRP(old.pop("NRP"))
    # alkaloid
    if "Alkaloid" in old:
        new["alkaloid"] = {}
        old_alk = old.pop("Alkaloid")
        if "alkaloid_subclass" in old_alk:
            new["alkaloid"]["subclass"] = old_alk.pop("alkaloid_subclass")
    # polyketide
    if "Polyketide" in old:
        pks = old.pop("Polyketide")
        if "Saccharide" in pks:  # BGC151
            old["Saccharide"] = pks.pop("Saccharide")
        new["polyketide"] = convert_pks(pks)
    # other
    if "Other" in old:
        old_other = old.pop("Other")
        sub = old_other.pop("other_subclass", old_other.pop("biosyn_class", None))
        if isinstance(sub, list):
            sub = sub[0]
        assert not old_other, old_other
        if sub in ["other", "None"]:
            sub = None
        if sub:
            new["other"] = {"subclass": sub}

    # ripp
    if "RiPP" in old:
        new["ripp"] = convert_ripp(old.pop("RiPP"))
    # saccharide
    if "Saccharide" in old:
        new["saccharide"] = convert_sacc(old.pop("Saccharide"))
    # terpene
    if "Terpene" in old:
        new["terpene"] = convert_terpene(old.pop("Terpene"))

    if "Other" in old:
        if old["Other"].get("other_subclass", "None") in ["None", "other"]:
            old.pop("Other")

    assert not old, old
    return new


def trim_evidence(data: JSON, mapping: Dict[str, str]) -> None:
    if "evidence" not in data:
        return
    invalid = ["Other", "N/A", "None"]
    evidence = data.pop("evidence")
    if isinstance(evidence, list):
        evidence = [mapping.get(item.strip(), item.strip()) for item in evidence if item not in invalid]
    else:
        evidence = mapping.get(evidence.strip(), evidence.strip())
    if evidence:
        data["evidence"] = evidence


def convert_ripp(old: JSON) -> JSON:
    def convert_precursor(old_pre: JSON) -> JSON:
        if not old_pre.get("gene_id"):
            return None
        new = {
            "gene_id": old_pre.pop("gene_id"),
        }
        rename_optionals([
            ("core_pept_aa", "core_sequence"),
            ("crosslinks", "crosslinks"),
            ("cleavage_recogn_site", "cleavage_recogn_site"),
            ("recogn_motif", "recognition_motif"),
        ], old_pre, new)
        commas_to_list(new, "cleavage_recogn_site")
        if isinstance(new.get("gene_id"), list) and new["gene_id"]:
            new["gene_id"] = new["gene_id"][0]
        old_pre.pop("foll_pept_len", None)
        old_pre.pop("lead_pept_len", None)
        old_pre.pop("peptidase", None)
        commas_to_list(new, "core_sequence")
        new_links = []
        for crosslink in new.get("crosslinks", []):
            rename_optionals([
                ("AA_pos_1", "first_AA"),
                ("AA_pos_2", "second_AA"),
            ], crosslink, crosslink)
            new_links.append(crosslink)
        assert not old_pre, old_pre
        return new
    new = {}
    rename_optionals([
        ("ripp_subclass", "subclass"),
        ("lead_pept_len", "lead_pept_len"),
        ("peptidase", "peptidase"),
    ], old, new)
    cyc = old.pop("lin_cycl_ripp", "").lower()
    if cyc == "cyclic":
        new["cyclic"] = True
    elif cyc == "linear":
        new["cyclic"] = False
    else:
        assert cyc == "", cyc

    precursors = []
    for pre in old.pop("precursor_loci", []):
        new_pre = convert_precursor(pre)
        if new_pre:
            precursors.append(new_pre)
    if precursors:
        new["precursor_genes"] = precursors

    assert not old, old
    return new


def convert_terpene(old: JSON) -> JSON:
    new = {}
    rename_optionals([
        ("terpene_subclass", "structural_subclass"),
        ("prenyl_transf", "prenyltransferases"),
        ("terpene_precursor", "terpene_precursor"),
        ("terpene_synth_cycl", "terpene_synth_cycl"),
        ("terpene_c_len", "carbon_count_subclass"),
    ], old, new)
    assert not old, old
    return new


def convert_sacc(old: JSON) -> JSON:
    def convert_transferase(data: JSON) -> JSON:
        if "gt_gene" not in data:
            print("discarding incomplete glycosyltransferase")
            return None
        new = {
            "gene_id": data.pop("gt_gene"),
            "specificity": data.pop("gt_specificity"),
        }
        if new["specificity"] == "None":
            new["specificity"] = "Unknown"
        if new["specificity"] == "Unknown":
            if "evidence_gt_spec" not in data:
                data["evidence_gt_spec"] = "Sequence-based prediction"  # TODO yuck
        else:
            if new["specificity"] == "Other":
                new["specificity"] = data.pop("other_gt_spec")
        evidence = data.pop("evidence_gt_spec", None)
        if evidence is None:
            evidence = ["Sequence-based prediction"]  # TODO more yuck
        if isinstance(evidence, str):
            evidence = [evidence]
        conversion = {
            "structure-based inference": "Structure-based inference",
        }
        new["evidence"] = evidence
        trim_evidence(new, conversion)

        assert not data, data
        return new

    new = {}

    assert not ("saccharide_subclass" in old and "Sugar_subclass" in old)
    rename_optionals([
        ("saccharide_subclass", "subclass"),
        ("Sugar_subclass", "subclass"),
    ], old, new)
    commas_to_list(new, "sugar_subclusters")
    transferases = []
    sugars = []
    for transferase in old.pop("gt_genes", []):
        old_sugars = transferase.pop("sugar_subcluster", "")
        assert isinstance(old_sugars, str), old_sugars
        sugars.append([sug.strip() for sug in old_sugars.strip().split(",")])
        trans = convert_transferase(transferase)
        if trans:
            transferases.append(trans)
    if transferases:
        new["glycosyltransferases"] = transferases
    if sugars:
        new["sugar_subclusters"] = sugars

    assert not old, old
    return new


def convert_moiety(old: JSON) -> JSON:
    new = {}
    subs = old.pop("moiety_subcluster", old.pop("subcluster", None))
    if subs and subs != "unknown":
        new["subcluster"] = subs
    name = old.pop("chem_moiety", "").lower()
    if name:
        key = f"{name}_chem_moiety"
        if key in old:
            new["moiety"] = old.pop(key)
        else:
            new["moiety"] = name
    if "moiety" not in new:
        return None  # BGC361
    assert not old, old
    return new


def convert_compounds(old: JSON) -> JSON:
    def convert_compound(old_compound: JSON) -> JSON:
        new = {
            "compound": old_compound.pop("compound"),
        }
        rename_optionals([
            ("chem_synonyms", "chem_synonyms"),
            ("mol_mass", "mol_mass"),
            ("mass_ion_type", "mass_spec_ion_type"),
            ("molecular_formula", "molecular_formula"),
            ("chem_struct", "chem_struct"),
            ("evidence", "evidence"),
            ("chem_act", "chem_acts"),
        ], old_compound, new)

        if "evidence_struct" in old_compound:
            if "evidence" not in new:
                new["evidence"] = []
            new["evidence"].extend(old_compound.pop("evidence_struct"))

        mol_mass = new.get("mol_mass")
        if mol_mass is not None and float(mol_mass) < 0:
            new.pop("mol_mass")
            mol_mass = None
        if mol_mass is not None and not isinstance(mol_mass, float):
            new["mol_mass"] = float(mol_mass)

        if new.get("mass_spec_ion_type") in ["None"]:
            new.pop("mass_spec_ion_type")

        if "chem_acts" in new:
            extras = old_compound.pop("other_chem_act", [])
            if isinstance(extras, str):
                extras = []
            new["chem_acts"].extend(extras)
        trimmed_acts = []
        for act in new.pop("chem_acts", []):
            if act.lower() == "unknown":
                continue
            trimmed_acts.append(act)
        if trimmed_acts:
            new["chem_acts"] = trimmed_acts

        targets = old_compound.pop("chem_target", [])
        if targets:
            if isinstance(targets, str):
                targets = [t.strip() for t in targets.split(",") if t != "unknown"]
            if targets:
                new["chem_targets"] = [{"target": target} for target in targets]

        moieties = old_compound.pop("chem_moieties", None)
        if moieties:
            converted = []
            for moiety in moieties:
                new_moiety = convert_moiety(moiety)
                if new_moiety:
                    converted.append(new_moiety)
            if converted:
                new["chem_moieties"] = converted

        dbs = []
        for db in old_compound.pop("databases_deposited", []):
            name = db.lower()
            old_id = old_compound.pop(f"{name}_id")
            if isinstance(old_id, str):
                old_id = old_id.strip()
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


def convert_genes(old: JSON) -> Tuple[JSON, List[Any]]:
    def convert_gene(old_gene: JSON) -> JSON:
        new = {}
        rename_optionals([
            ("gene_name", "name"),
            ("mut_pheno", "mut_pheno"),
            ("gene_id", "id"),
            ("gene_annotation", "product"),
            ("tailoring", "tailoring"),
            ("gene_comments", "comments"),
        ], old_gene, new)
        name = new.get("name")
        if name is not None and name in ["", "No gene ID"]:
            new.pop("name")
        commas_to_list(new, "tailoring")
        if "gene_function" in old_gene:
            new_function = {}
            rename_optionals([
                ("gene_function", "category"),
                ("evidence_genefunction", "evidence"),
            ], old_gene, new_function)
            if not new.get("category", None):
                # don't go on to check evidence for a function that doesn't exist
                return new
            evidence = new_function.pop("evidence", [])
            if new_function["category"] == "Unknown":
                evidence = []  # no function? no evidence!
            evidence = [item.strip() for item in evidence if item not in ["Other", "N/A"]]
            if evidence:
                new_function["evidence"] = evidence
                new["functions"] = [new_function]
            else:
                print("discarding gene function from", new.get("name", new["id"]), "due to missing evidence")

        external = old_gene.pop("not_in_gbk", False)
        if not external:
            # discard any positions, since the genbank should be the real thing
            old_gene.pop("gene_startpos", None)
            old_gene.pop("gene_endpos", None)
        else:
            start = old_gene.pop("gene_startpos", None)
            end = old_gene.pop("gene_endpos", None)
            if start is None or end is None:
                warnings.warn(f"discarding external gene with no location: {new['name']}")
                return None
            warnings.warn(f"unhandled extra gene: {new.get('name', gene.get('gene_id', '???'))}")  # TODO
            return None

        assert not old_gene, old_gene
        return new

    new = {}
    publications = []
    genes = []
    for gene in old.pop("gene", []):
        if "gene_pubs" in gene:
            pubs = gene.pop("gene_pubs")
            for pub in pubs:
                if not pub:
                    continue
                if pub == "None":
                    continue
                if pub.isdigit():
                    publications.append(f"pubmed:{pub}")
                    continue
                assert False, pub
        new_gene = convert_gene(gene)
        if new_gene is not None:
            genes.append(new_gene)
    if genes:
        new["annotations"] = genes
    if "operon" in old:
        operons = old.pop("operon")
        valid = []
        for operon in operons:
            op_genes = operon.pop("operon_genes")
            if op_genes:
                operon["genes"] = op_genes
            rename_optionals([
                ("evidence_operon", "evidence"),
            ], operon, operon)
            commas_to_list(operon, "evidence")
            evidence = [ev for ev in operon.pop("evidence", []) if ev != "Other"]
            if not evidence:
                print("discarding operon", operon, "due to missing evidence")
                continue
            operon["evidence"] = evidence
        new["operons"] = [op for op in valid if op]

    assert not old, old
    return new, publications


def rename_optionals(keys: List[Tuple[str, str]], old: JSON, new: JSON) -> None:
    for src, dst in keys:
        if src in old:
            val = old.pop(src)
            if val:
                new[dst] = val


def convert_NRP(old: JSON) -> JSON:
    def convert_module(old_module: JSON) -> JSON:
        new = {}
        rename_optionals([
            ("module_nr", "module_number"),
            ("cdom_subtype", "c_dom_subtype"),
        ], old_module, new)

        if new.get("module_number").lower() == "x":
            new.pop("module_number")

        if new.get("c_dom_subtype") in ["N/A", "None"]:
            new.pop("c_dom_subtype")
        if new.get("c_dom_subtype") in ["Other"]:
            new["c_dom_subtype"] = "Unknown"

        assert not old, old
        return new

    new = {}
    rename_optionals([
        ("lipid_moiety", "lipid_moiety"),
        ("subclass", "subclass"),
        ("nrps_release_type", "release_type"),
    ], old, new)

    commas_to_list(new, "release_type")

    if "lin_cycl_nrp" in old:
        assert old["lin_cycl_nrp"] in VALID_CYCLICS
        new["cyclic"] = old.pop("lin_cycl_nrp") == "Cyclic"

    thios = []
    if "nrps_thioesterase" in old or "nrps_te_type" in old:
        thio_type = old.pop("nrps_te_type")
        if thio_type in ["Both", "Other", "other", "None"]:
            thio_type = "Unknown"
        for thio in old.pop("nrps_thioesterase", []):
            thios.append({
                "gene": thio,
                "thioesterase_type": thio_type,
            })
    if thios:
        new["thioesterases"] = thios

    genes = []
    for gene in old.pop("nrps_genes", []):
        if not gene:
            continue
        modules = gene.pop("modules", [])
        if not modules and "nrps_module" in gene:
            modules = gene.pop("nrps_module")
        if modules:
            modules = [convert_module(module) for module in modules]
            gene["modules"] = sorted(modules, key=lambda module: module.get("module_number", "ZZZ"))
        rename_optionals([
            ("nrps_gene", "gene_id"),
        ], gene, gene)
        genes.append(gene)
    if genes:
        new["nrps_genes"] = genes

    assert not old, old
    return new


def convert_pks_synthase(old: JSON) -> JSON:
    def convert_module(old_module: JSON) -> JSON:
        new = {}
        rename_optionals([
            ("module_nr", "module_number"),
            ("gene", "genes"),
            ("pks_domains", "domains"),
            ("at_substr_spec", "at_specificities"),
            ("evidence_at_spec", "evidence"),
            ("pks_mod_doms", "pks_mod_doms"),
            ("kr_stereochem", "kr_stereochem"),
        ], old_module, new)
        commas_to_list(new, "genes")
        commas_to_list(new, "at_specificities")

        if new.get("evidence") == "Other":
            new.pop("evidence")

        stereo = new.get("kr_stereochem")
        if stereo == "B-group":
            new["kr_stereochem"] = "L-OH"   # TODO this is a bit iffy feeling
        elif stereo == "A-group":
            new["kr_stereochem"] = "D-OH"   # TODO this is a bit iffy feeling

        # remove modification domains if the old value was just saying there were none
        mods = new.get("pks_mod_doms")
        if not mods or mods == "None":
            new.pop("pks_mod_doms")

        non_canonical = old_module.pop("pks_mod_skip_iter", "Neither")
        non_canonical_evidence = old_module.pop("pks_evidence_skip_iter", None)
        mapping = {
            "Iterated": "iterated",
            "Skipped": "skipped",
            "Non-elongating": "non_elongating",
        }
        if non_canonical != "Neither":
            new["non_canonical"] = {}
            new["non_canonical"][mapping[non_canonical]] = True
        if non_canonical_evidence is not None and non_canonical_evidence not in ["Other", "None"]:
            if new.get("evidence") is None:
                new["evidence"] = non_canonical_evidence  # TODO also iffy
            elif non_canonical_evidence != new["evidence"]:
                new["non_canonical"]["evidence"] = non_canonical_evidence
                commas_to_list(new["non_canonical"], "evidence")
        assert not old_module, old_module
        return new

    new = {}
    new["genes"] = []
    new["genes"].extend(gene["mod_pks_gene"] for gene in old.get("mod_pks_genes", []))
    new["genes"].extend(old.pop("pks_genes", []))
    modules = []
    for gene in old.pop("mod_pks_genes", []):
        old_modules = gene.pop("pks_module", [])
        name = gene.pop("mod_pks_gene")
        for module in old_modules:
            module["gene"] = name
            modules.append(convert_module(module))
    if modules:
        new["modules"] = sorted(modules, key=lambda module: module["module_number"])

    rename_optionals([
        ("pufa_mod_doms", "pufa_modification_domains"),
        ("pk_subclass", "subclass"),
    ], old, new)
    commas_to_list(new, "subclass")
    thios = []
    thio_type = old.pop("pks_te_type")
    if thio_type in ["Both", "Other"]:
        thio_type = "Unknown"
    for thio in old.pop("pks_thioesterase", []):
        thios.append({
            "gene": thio,
            "thioesterase_type": thio_type,
        })
    if thios:
        new["thioesterases"] = thios

    if any(key in old for key in ["nr_iterations", "iterative_subtype", "iter_cycl_type"]):
        iterative = {
            "nr_iterations": old.pop("nr_iterations"),
            "subtype": old.pop("iterative_subtype"),
            "cyclization_type": old.pop("iter_cycl_type", "Unknown"),
        }
        new["iterative"] = iterative

    trans_at = old.pop("trans_at", None)
    if trans_at:
        new["trans_at"] = {"genes": trans_at}

    assert not old, old
    return new


def commas_to_list(data: JSON, key: str) -> None:
    if key not in data:
        return
    value = data[key]
    if isinstance(value, list):
        return
    data[key] = [val.strip() for val in value.split(",")]


def convert_pks(old: JSON) -> JSON:
    new = {}
    rename_optionals([
        ("starter_unit", "starter_unit"),
        ("pks_subclass", "subclasses"),
        ("pks_release_type", "release_type"),
    ], old, new)

    # merge pk_subclass and pks_subclass
    if "pk_subclass" in old:
        new["subclasses"].append(old.pop("pk_subclass"))

    commas_to_list(new, "release_type")
    if "lin_cycl_pk" in old:
        assert old["lin_cycl_pk"] in VALID_CYCLICS
        new["cyclic"] = old.pop("lin_cycl_pk") == "Cyclic"

    if "mod_pks_genes" in old or "pks_genes" in old:  # genes are required in new synthase annotations
        synth = convert_pks_synthase(old)
        if synth["genes"]:
            new["synthases"] = [synth]

    thio = old.pop("pks_te_type", None)
    if thio and thio != "None":
        assert False, thio

#       self.cyclases = raw.get("cyclases")  # list[str]
#        self.ketide_length = raw.get("ketide_length")  # int
    assert not old, old
    return new


def transform(data: JSON, mibig_version: str) -> JSON:
    transformed = {}
    transformed["cluster"] = build_cluster(data.pop("general_params"))

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
    warnings.warn("still using dummy changelog")  # TODO
    data.pop("changelog", None)
    transformed["changelog"] = changelog

    # discard embargo
    data.pop("embargo", None)

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

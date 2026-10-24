# Draft maintainer note: historical GtoP namespace scope and 11 current identity conflicts

Prepared 24 October 2026. **Not sent or posted.** All evidence is public; no patient data is involved.

Thank you for fixing the GtoP target-ID substitution in [PR624](https://github.com/dgidb/dgidb-v5/pull/624). We independently measured its historical scope and verified that this exact pattern is absent from the tested June and October 2026 exports. This note is about downstream cleanup and a separate set of residual identity conflicts.

**Reproducible historical scope**

In the official DGIdb **2024-12** interactions export, 13,146 distinct gene–GtoP-ligand pairs meet a strict numeric-swap criterion. These are 12,994 gene–normalized-drug pairs across 1,403 genes, out of 16,365 GtoP rows. The criterion requires the assigned HGNC gene to equal the human NCBI mapping of the mistakenly substituted integer and to disagree with the GtoP target's HGNC identity. It excludes 1,158 other identity conflicts. A second implementation reproduces the result.

The original bad line was:

```ruby
create_gene_claim("NCBIGENE:#{line['Target ID']}", GeneNomenclature::NCBI_ID)
```

The existing repair correctly uses a GtoP claim dictionary instead. Named regression cases are `701` (CLCNKB versus BUB1B), `2760` (PDCD1 versus GM2A), `2628` (MS4A1 versus GATM), `215` (DRD2 versus ABCD1), `16` (CHRM4 versus AARS1), `24` (ADRA1D versus ABCA4) and `1376` (PTGS2 versus CPT2). [Pair-level census](../results/collision_audit/identity_conflicts.csv), [named evidence](../results/collision_audit/named_cases.csv).

**Current residuals requiring investigation**

| GtoP target | Source human identity | Assigned normalized concept | GtoP edges in October export and live API |
|---|---|---|---:|
| 1310 | PDE10A; HGNC:8772; Entrez:10846 | LINC00473; hgnc:21160 | 8 |
| 2847 | MCL1; HGNC:6943; Entrez:4170 | ADAMTSL4-AS1; hgnc:32041 | 3 |

Their current claim namespace is already `IUPHAR.RECEPTOR`. They are **not classified as proven instances of the old numeric substitution**. We have not established whether the cause is stale aliases, a resolver merge or another path. The live API returns null for the intended `hgnc:8772` and `hgnc:6943` concept queries. We did not adjudicate every other-source interaction returned for the incorrect concepts.

Read-only reproduction at `https://dgidb.org/api/graphql`:

```graphql
{
  assignedPDE10A: gene(conceptId: "hgnc:21160") {
    name
    interactions {
      sources { fullName }
      drug { name conceptId }
      interactionTypes { type }
    }
  }
  assignedMCL1: gene(conceptId: "hgnc:32041") {
    name
    interactions {
      sources { fullName }
      drug { name conceptId }
      interactionTypes { type }
    }
  }
  expectedPDE10A: gene(conceptId: "hgnc:8772") { name }
  expectedMCL1: gene(conceptId: "hgnc:6943") { name }
}
```

Filter interactions to `sources.fullName == "Guide to Pharmacology"`. The eight PDE10A ligands are mardepodect, PDE4 inhibitor 16, PBF-999, compound 19 [PMID:22142545], TC-E 5005, compound 27 [PMID:22142545], papaverine and MK-8189. The three MCL1 ligands are venetoclax, AZD5991 and obatoclax. Identity assignment is the finding; this note does not certify affinity or selectivity. [Complete export rows](../results/collision_audit/current_identity_conflicts.csv), [original live query, response evidence and hashes](../results/collision_audit/live_api_evidence.json).

**Suggested repair and prevention**

Keep the PR624 importer repair. Add a validation gate **after normalization**, because a correctly namespaced source claim can still end up in the wrong normalized concept. Reference logic:

```python
source_id = ("IUPHAR.RECEPTOR", row.target_id)
source_target = gtop_targets[source_id]
assert source_target.species == "Homo sapiens"

expected = validate_human_gene_cross_references(source_target)
resolved = normalize(source_id)
if not expected or not resolved or resolved.hgnc not in expected.hgnc_ids:
    quarantine(row, source_target, resolved, reason="target_identity_conflict")
else:
    persist_with_provenance(row, source_target, resolved)
```

This is reference logic, not a claimed drop-in patch for the Rails application. A production implementation must handle protein complexes and multi-gene targets explicitly. Do not resolve a disagreement by choosing whichever symbol happens to match. Namespace parsing must also reject foreign IDs: the current GtoP target 3348 HGNC field contains `MGI:103220|7013`. No audited DGIdb row uses this target, so it contributes zero to our counts; the parser quarantines the MGI token rather than manufacturing an HGNC identifier.

For the residual clusters, trace alias provenance and rebuild affected concepts from clean assertions. Do not copy an entire mixed concept onto a new gene or merely rename it. Re-export and invalidate downstream caches after verifying each source assertion.

Acceptance checks:

1. Every imported human source assertion reconciles with its target cross-references or is explicitly quarantined.
2. The seven historical sentinel collisions cannot resolve to their NCBI lookalikes.
3. The 11 listed residuals resolve to the source genes or leave the published assertion set pending review.
4. Known correct target–ligand pairs survive; deletion of all source rows is not a passing repair.
5. A release manifest exposes source versions, correction date, row counts and an invalidated-edge list.

A correction notice for historical exports would help consumers. A 2025 paper's Table 1 still reports DGIdb-derived GM2A interactions with nivolumab, cemiplimab and pembrolizumab. [PMID 41165052; DOI 10.1002/epi4.70175](https://pmc.ncbi.nlm.nih.gov/articles/PMC12716302/). This establishes persistence in a derivative; it does not measure clinical use or harm.

**Audit package**

Run `python3 scripts/replay_collision_audit.py` with the preserved public cache. It verifies 52 inputs and reproduces eight output files byte for byte. The six offline behavioral tests run with `uv run pytest tests/test_collision_audit.py`. Source URLs, request bodies, retrieval times and SHA-256 hashes are in the [manifest](../results/collision_audit/input_manifest.json). Current GtoP/NCBI reference identities are used; we did not replay the historical Rails database or recover the complete 2024 upstream reference snapshot.

The result is source-specific. It does not establish that GtoP's primary pharmacology assertions are wrong or that all DGIdb sources share this defect. Our separate complete ChEMBL mechanism-collection comparison found no overlap with the 7,800 mappable historical wrong gene–compound pairs.

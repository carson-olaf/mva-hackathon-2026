# Burrell supplementary-methods check

24 October 2026. Public-source check within the root's 20-minute cap (started 20:32:21 UTC; hard stop 20:52:21 UTC). Completed early after retrieving and examining the publisher's supplementary PDF. No patient/private files, contacts, raw-source files saved to disk or Git operations. HTTPS response bytes and PDF text remained in memory; this note is the only file written by this agent.

**Result: the supplement was retrieved, but it does not resolve uridine-alone sufficiency or the dose assignment.** The 45–55% result must remain attributed to nucleoside supplementation as defined by the four-component main methods, not to uridine alone. The supplement also does not prove that every component is necessary. No component-deletion or uridine-only arm was found.

## Primary source obtained

- Article: Burrell et al., PMID **23446422**, DOI **10.1038/nature11935**.
- Publisher article URL: https://www.nature.com/articles/nature11935
- Actual supplementary PDF href, extracted from a successful HTTP 200 article response:
  https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11935/MediaObjects/41586_2013_BFnature11935_MOESM68_ESM.pdf
- PDF response: HTTP 200; application/pdf; **5,145,921 bytes**.
- SHA-256 of those ephemeral response bytes: **8bceb7b460d597c3ad0c667973a56c8de942c0b36dc5844439acc18258934ed6**.
- Text extraction: installed `/opt/homebrew/bin/pdftotext -layout - -`, passing PDF bytes on stdin and reading stdout; no source file was created.
- **42 nonempty PDF pages**. Figures/tables/reference material occupy pages 1–16; a Sweave analysis document begins on page 17.

## What the retrieved supplement establishes

1. **Supplementary Figure 13, PDF page 13.** Panels compare untreated cells with nucleoside supplementation. The legend specifies **48 hours**. Panel a measures prometaphase DNA damage in HT29, HT55, SW620 and SW1116; panels b/c examine segregation errors in CIN-negative lines and NCIH508. It supplies neither a concentration nor component-specific treatment arms.
2. **Supplementary Figure 14, PDF page 14.** The treatment labels remain untreated versus nucleosides. These panels check proliferation, cell-cycle distribution, mitotic index and ATP. The legend supplies experimental duration/controls but not a nucleoside dose or uridine-only result.
3. A full-file text search found nucleoside mentions **only on pages 13 and 14**, and **no occurrence of uridine**. Numeric occurrences of 0.3 elsewhere are not evidence of a treatment concentration: the relevant Figure 14 occurrence is a proliferation-rate axis, and later occurrences belong to analysis code/plots. No dose can be inferred from those numbers.
4. The Sweave material concerns genomic/copy-number statistical analyses; its existence does not supply a hidden single-nucleoside treatment result.

**Still unverified after obtaining the supplement:** the exact allocation of the main-method concentrations **0.3 or 30 µM** to the Figure 4 rescue experiments; whether the intended notation is explicitly per component; an effective concentration/duration for **uridine alone**. The main paper's A/C/G/U preparation remains the verified mixture definition from the prior primary-text inspection. Do not substitute a newly retrieved document for a result it never reports.

The decisive remaining experiment is therefore a direct uridine-only versus full-mixture comparison, followed by exact-genotype rescue and harm testing. A mixture-only published comparison cannot establish either single-component sufficiency or four-component necessity.

## Correction checked

- Primary publisher correction: https://www.nature.com/articles/nature12369 ; PDF https://www.nature.com/articles/nature12369.pdf ; DOI **10.1038/nature12369**.
- The indexed primary Nature correction text was read. It adds omitted **TCGA acknowledgements and dbGaP accession information**. It contains no change to nucleoside composition, concentration or rescue endpoints.
- Direct web-tool opens of both correction URLs failed through the Nature identity redirect. Thus the correction's indexed primary full text, rather than its rendered PDF, is the verified access scope.

## Retrieval limits and failed routes

- The web tool failed to open the Nature article/correction via its identity redirect (HTTP 500), the true media.springernature.com supplement URL (inaccessible), legacy supplementary URLs (cache miss/inaccessible), and the UCL thesis PDF (cache miss).
- An in-memory request with the repository's TLS-capable `requests` client succeeded. The actual supplement is **MOESM68**, not the previously guessed MOESM1; the legacy `nature11935-s1.pdf` returned HTTP 404.
- UCL's article repository record exposes the main manuscript, not an independently useful supplementary file. The public doctoral-thesis repository record (https://discovery.ucl.ac.uk/id/eprint/1383596/) was read, but its PDF was not obtained and no thesis-derived dose or efficacy claim is made.
- No statement is made about unpublished laboratory records or unavailable raw experimental treatment sheets. No author or family contact occurred.

## Minimal public reproduction

Use a current TLS-capable Python environment with requests and Poppler pdftotext. Fetch the article HTML in memory, extract its PDF href containing `MOESM68`, fetch that response in memory, verify the hash above, and pass bytes to `pdftotext -layout - -`. Examine PDF pages 13–14 and search the entire text for nucleoside, uridine and the candidate doses. The saved report records observed source identity and scope; it does not claim that raw bytes were archived.

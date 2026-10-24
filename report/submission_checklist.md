# Submission checklist — carcango

Prepared 24 October 2026. **Deadline: 23:59 UTC / 16:59 Vancouver today.** This is a local handoff, not a submission receipt. No push, upload, repository publication or submission was performed. [Official timeline](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py).

**The proposal:** approved uridine triacetate as a candidate for laboratory testing of a nucleotide-responsive component of chromosome instability. The first page states the candidate, two unproved biological transfers and the rejection test. The published mixture result remains distinct from unproved uridine-alone efficacy; the same-range DNA-damage hazard stays prominent.

## Ready locally

- [x] [Track 2 report PDF](carcango_track2_report.pdf), with the required username in its filename. Editable source: [track2_report.md](track2_report.md).
- [x] [One-page evidence summary PDF](../results/carcango_codex_summary.pdf), with [Markdown source](../results/codex-summary.md). This is a companion; the form does not require a separate summary upload.
- [x] [Three-minute pitch script](pitch_script.md). Script only; actual recorded duration and audio quality remain unverified.
- [x] [Track 1 CSV](../results/track1_carcango_hpo_recessive.csv): 12 official columns, GRCh38, chr-prefixed chromosomes, one primary pair. The file uses the VCF sample name, WGS_EX2312012.
- [x] [Track 1 methods PDF](carcango_track1_methods.pdf), with username-bearing filename.
- [x] Required acknowledgement retained verbatim; OpenAI / Pro / training and data-sharing disabled disclosure retained. Settings are reported by Carson, not inferred from account inspection.
- [x] Burrell publisher supplement obtained; [audit](../results/burrell_supplement_audit.md) records the hash, access route and remaining dose/component questions.
- [x] All eight external-review corrections applied in separate local commits. Pandey and both Burrell transfer limitations verified; AI disclosure and chronology preserved.
- [x] Reproducible exposure arithmetic and public source audit included. No uridine-only efficacy, continuous exposure, cancer prevention or clinician signoff is claimed.

Final file hashes, page counts, test result and local checks are recorded in [submission validation](../results/submission_validation.json). This is separate from the earlier computational-verification snapshot; that older manifest is not relabeled as a current run.

## Carson's actions before the deadline

1. **Read the report opening and rehearse the pitch.** Record a video no longer than three minutes. Upload it to **YouTube or Vimeo** and verify the link is accessible to judges. No video was created or checked in this run.
2. **Provide the GitHub repository URL.** `git remote -v` returned no configured remote in this checkout. No remote URL, publication or judge access has been verified. The code is prepared locally; creating/pushing the remote remains a separately authorized action. Do not upload `data/`, `private/`, `work/`, raw reads, alignments or genotype-scale tables.
3. **Confirm the registered account and team submitter.** Use the intended Hugging Face account. The live account's quota and previous receipts were not inspected. Track 2 allows at most **three** submissions and judges the **latest**; do not replace a final entry with an older draft.
4. **Submit Track 2:** upload `carcango_track2_report.pdf`, enter the GitHub URL and the YouTube/Vimeo URL, and check the report title/team. Retain the actual confirmation. The PDF includes figures and external citations; local companion paths refer to the separately supplied repository.
5. **Submit Track 1 if claiming that entry:** use `track1_carcango_hpo_recessive.csv`, `carcango_track1_methods.pdf` and the GitHub URL. At most **six** submissions are allowed; the highest score is shown. Production acceptance of the proband ID is unverified. Do not replace it with `PROBAND01` solely from the public code's error string: that value is the local-development fallback, while production loads a private mapping. Preserve the provenance-backed ID unless organizer evidence establishes another alias.
6. **Check both receipts and access.** Confirm the uploaded report is the final PDF, the video plays, and the repository location is correct. A local file, passing test or completed form is not proof of submission.

Official sources: [Track 2 form and instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/submit_track2.py), [Track 1 form](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py), [public ground-truth loader](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/groundtruth.py), [FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py). Only public code was inspected; no private answer key was accessed.

## After close — separate human-authorized work

- The repository may be private during the event but must become public once the hackathon ends. Review files and history before authorizing release. No push or publication is authorized by this checklist.
- Preserve source-specific licenses. Original submission materials use CC BY 4.0; restricted patient data and separately licensed third-party material do not inherit that license.
- The exact dataset citation required at publication remains unresolved after inspecting the [public Synapse page](https://www.synapse.org/Synapse:syn76251147/wiki/642892). Obtain it before publication; do not invent a DOI.
- Follow the manuscript embargo until the organizer summary/preprint. Conference preliminary findings require prior written approval.
- Official Rules require deletion of underlying data, intermediates and derived datasets within **30 days of close**, followed by confirmation to **cnb@sagebase.org**. Arrange a rules-aware retention/deletion review; this task did not delete data or send email.

[Official Rules, data use, citation, acknowledgement and publication](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py). No raw patient data should be packaged with this handoff.

# Reporter CDS source audit

**Audit date:** 2026-09-30
**Status (updated 2026-10-06):** Exact pJ251-GERC sequence comparison is not completed and will be skipped. The user reports that Addgene registration requires organization information they cannot provide after graduation; the user has decided this is an objective access constraint. Keep the CDS reference provisional and disclose this limitation.

## Evidence checked

- The current [provisional 90-nt FASTA](../../data/reference/sfgfp_first90nt.fasta) is attributed to Appendix C.1.1 of Gilliot's thesis. Its sequence starts with `ATG` and supplies 87 coding nucleotides after the start codon for the ViennaRNA context.
- The current [Addgene sequence page](https://www.addgene.org/47441/sequences/) lists one depositor-provided full sequence plus three Addgene-verified partial sequences and requires login to download or copy them. Addgene's [current help page](https://help.addgene.org/hc/en-us/articles/44210206209549-Are-there-any-requirements-for-viewing-plasmid-sequence-data-on-the-Addgene-website) says accounts are free and sequence files are available from the “Sequence” link after login; some records may require acknowledging the Affinity Reagent Sequence Policy. The record's full sequence is depositor-provided and may differ from Addgene sequencing results, so provenance must be recorded when retrieved.
- Kosuri et al. identify pGERC as the sfGFP reporter construct and report deposition of the sequence to Addgene, but the accessible paper text does not supply the full reporter ORF in a form that verifies these exact 90 nt.

## Sequence checksum

SHA-256 of the current FASTA file: `684f537fad73723f43efe393665386b50dc6261105095396e4569469510865cc`.

## Conclusion and next action

The current 90-nt sequence remains suitable for method development but is **not verified against the assayed pJ251-GERC construct**. Do not describe it as the exact construct sequence. The Addgene route is skipped by user decision because registration is inaccessible under the user's current organization circumstances. Preserve the provisional label and report the unresolved sequence provenance as a limitation. If an authorized sequence file becomes available through another route later, extract the sfGFP ORF in the correct orientation and compare the first 90 nt, start codon, and cloning junction while recording its provenance.

References: [Addgene pJ251-GERC sequence page](https://www.addgene.org/47441/sequences/); [Kosuri et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3752251/); [Gilliot thesis](https://research-information.bris.ac.uk/files/432056523/Final_Copy_2024_01_23_Gilliot_P_PhD_Redacted.pdf).

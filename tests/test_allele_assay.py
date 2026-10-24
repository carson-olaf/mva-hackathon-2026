from mva_hackathon_2026.allele_assay import digest, site_peptides


def test_new_lysine_creates_cut_and_one_missed_product_but_proline_blocks_it():
    # A synthetic protein exposes the assay's central mutant-cleavage confound.
    wt, mutant = 'AAANAAKPEER', 'AAAKAAKPEER'
    assert [p['sequence'] for p in site_peptides(wt, 4, max_missed=0)] == [wt]
    products = site_peptides(mutant, 4)
    assert [(p['sequence'], p['missed_cleavages']) for p in products] == [
        ('AAAK', 0), (mutant, 1)]
    assert [p['sequence'] for p in digest(mutant, max_missed=0)] == ['AAAK', 'AAKPEER']
    assert [p['sequence'] for p in site_peptides(mutant, 4, cleavage_after='E',
                                               block_before_proline=False,
                                               max_missed=0)] == ['AAAKAAKPE']

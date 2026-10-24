from mva_hackathon_2026.track1 import annotate, nmd, phenotype_prior, rank_models


def test_gene_blind_recessive_pair_beats_single_missense():
    text='HP:1\tRare\t1\tGENEA\tdisease\nHP:2\tBroad\t1\tGENEA\tdisease\nHP:2\tBroad\t2\tGENEB\tdisease'
    prior, _, _=phenotype_prior(text,['HP:1','HP:2'])
    calls=[{'gene':'GENEA','pos':1,'het':True,'consequence':'stop_gained'}, {'gene':'GENEA','pos':10,'het':True,'consequence':'missense_variant'}, {'gene':'GENEB','pos':20,'het':True,'consequence':'missense_variant'}]
    ranked=rank_models(calls,prior)
    assert ranked[0]['gene']=='GENEA'
    assert ranked[0]['model']=='recessive_compound_het_unphased'
    calls[0]['max_population_af']=0.2
    assert all(len(r['variants'])==1 for r in rank_models(calls,prior))


def test_coding_reference_strand_and_nmd_boundary():
    # Two real coding rules tested with a synthetic transcript, not patient data.
    model={'positions':list(range(101,110)),'sequence':'ATGAAATAA','strand':1,'exon_numbers':[1]*9,'exons':[{'start':101,'end':109},{'start':201,'end':209}],'transcript':'SYNTHETIC.1'}
    assert annotate(model,104,'A','T')['consequence']=='stop_gained'
    assert annotate(model,104,'C','T')=={'error':'reference_mismatch'}
    assert not nmd(model,{'coding_position':4})['nmd_expected_50_55_rule']
    reverse={**model,'positions':list(range(109,100,-1)),'strand':-1}
    assert annotate(reverse,106,'T','A')['consequence']=='stop_gained'
    long_model={**model,'exons':[{'start':101,'end':180},{'start':201,'end':209}]}
    assert nmd(long_model,{'coding_position':4})['bases_after_stop_before_final_junction']==74
    assert nmd(long_model,{'coding_position':4})['nmd_expected_50_55_rule']

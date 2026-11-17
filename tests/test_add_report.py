"""Dry runs account for input rows, including partial resolution and expansion."""
import pandas as pd
from pharos_gap_kg import add, ingest


class Graph:
    ti = {'UniProtKB:P1': 0, 'UniProtKB:P2': 1}
    di = {'MONDO:1': 0}


def symbols(value):
    return {'A': ['UniProtKB:P1'], 'B': ['UniProtKB:P2'],
            'AMB': ['UniProtKB:P1', 'UniProtKB:P2']}.get(value, [])


def proposal():
    return {'map': True, 'subject': {'columns': ['gene'], 'resolver': 'gene_symbol',
                                    'separator': '/', 'category': 'biolink:Protein'},
            'object': {'constant': 'MONDO:1', 'category': 'biolink:Disease'},
            'predicate': 'biolink:genetically_associated_with', 'association_category': 'biolink:Association',
            'knowledge_level': 'statistical_association', 'agent_type': 'data_analysis_pipeline'}


def report(df, p, accepted=True):
    transformed, rows = add.transform(p, df, symbols, None, 'infores:test')
    v = {'accepted': accepted, '_p2': transformed, '_rows': rows,
         'rules': [{'checks': {}} for _ in ingest.rules_of(transformed)]}
    return add.mapping_report(v, df, p, Graph(), None, symbols, 2)


def test_partial_and_multiple_matches_are_visible_and_rows_count_once():
    df = pd.DataFrame({'gene': ['A/UNKNOWN', 'AMB', 'MISSING', None]})
    out = report(df, proposal())
    assert out['source_rows'] == out['accounted_source_rows'] == 4
    assert out['outcomes'] == {'graph_pair': 2, 'unresolved': 2}
    assert out['issue_counts']['unmatched_symbol'] == 2
    assert out['issue_counts']['multiple_symbol_matches'] == 1
    assert any(x['value'] == 'UNKNOWN' for x in out['issue_examples'])
    assert out['outcome_examples']['graph_pair'][0] == {'source_row': 1, 'sheet_row': 3}
    assert out['exported_edges'] is None


def test_gaps_skips_and_unmapped_split_values_are_separate():
    p = proposal()
    p.update(split_by='kind', by_value={'use': {}, 'omit': {'skip': 'not relevant'}})
    p['keep_if'] = {'column': 'p', 'at_least': 0.05, 'reason': 'weak'}
    df = pd.DataFrame({'gene': ['A', 'B', 'A', 'A'], 'kind': ['use', 'use', 'omit', 'other'],
                       'p': ['0.01', '0.5', '0.01', '0.01']})
    out = report(df, p, accepted=False)
    assert out['outcomes'] == {'gap': 1, 'graph_pair': 1, 'skipped': 1, 'unmapped_split_value': 1}
    assert not out['accepted']
    assert out['accounted_source_rows'] == 4


def test_refused_and_declined_proposals_do_not_promise_an_export():
    df = pd.DataFrame({'gene': ['A', 'B']})
    for v, p, outcome in [({'accepted': False}, {}, 'not_evaluated'),
                           ({'accepted': True, 'map': False}, {'map': False}, 'declined')]:
        out = add.mapping_report(v, df, p, Graph(), None, symbols, 1)
        assert out['outcomes'] == {outcome: 2}
        assert out['exported_edges'] is None


def test_issue_examples_are_bounded_but_counts_are_not():
    df = pd.DataFrame({'gene': ['UNKNOWN'] * 30})
    out = report(df, proposal())
    assert out['accounted_source_rows'] == 30
    assert out['issue_counts']['unmatched_symbol'] == 30
    assert len(out['issue_examples']) == 10
    assert out['issue_examples_truncated']

"""Full-table diagnostics are observational and must not alter imported data."""
import pandas as pd
from pharos_gap_kg import add, ingest


def test_profile_reports_missingness_and_repeated_rows_without_dropping_them():
    df = pd.DataFrame({'gene': ['APOE', 'APOE', None, 'NA'], 'value': ['1', '1', '2', '3']})
    before = df.copy(deep=True)
    out = add.table_profile(df)
    pd.testing.assert_frame_equal(df, before)
    assert out['rows'] == 4
    assert out['duplicate_rows_beyond_first'] == 1
    for col in out['columns']:
        expected = ingest.present(df[col['name']])
        assert col['rows_with_a_value'] == int(expected.sum())
        assert col['missing_rows'] == int((~expected).sum())
    assert out['columns'][1]['distinct_values'] == 3


def test_profile_counts_the_whole_table_not_only_preview_rows():
    df = pd.DataFrame({'gene': ['APOE'] * 20 + [None]})
    out = add.table_profile(df)
    assert out['columns'][0]['missing_rows'] == 1
    assert out['rows'] == 21
    assert out['duplicate_rows_beyond_first'] == 19


def test_selected_sheet_card_exposes_profile(tmp_path):
    class Graph:
        di = {}
        def find(self, *args):
            return {'matches': []}
    class Session:
        G = Graph()
    source = tmp_path / 'genes.tsv'
    source.write_text('gene\tvalue\nAPOE\t1\nAPOE\t1\n')
    out = add.tool_card(Session(), str(source), add.STUDIES['A']['about'], source.name, 1)
    assert out['profile']['duplicate_rows_beyond_first'] == 1
    assert out['stage'] == 'mapping'

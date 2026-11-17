"""Malformed MCP requests must fail before opening files or writing sidecars."""
import pytest
from pharos_gap_kg import add


@pytest.mark.parametrize('field,value', [('doi', ''), ('doi', '...'), ('citation', None),
                                         ('licence', 42), ('disease', []), ('published', '2026-02-30')])
def test_invalid_metadata_is_refused_before_file_access(field, value):
    about = {**add.STUDIES['A']['about'], field: value}
    for result in (add.tool_card(None, '/missing.xlsx', about),
                   add.tool_check(None, '/missing.xlsx', about, 'mapping', {}),
                   add.tool_join(None, '/missing.xlsx', about, 'sheet', 1, {})):
        assert result['problems']
        assert not result.get('accepted', result.get('joined', False))


def test_unknown_stage_is_not_treated_as_mapping():
    result = add.tool_check(None, '/missing.xlsx', add.STUDIES['A']['about'], 'typo', {})
    assert result == {'accepted': False, 'problems': ['stage must be sheet or mapping']}


def test_registered_study_descriptions_remain_accepted():
    for study in add.STUDIES.values():
        assert add._about_problems(study['about']) == []

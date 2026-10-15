"""The person test's frozen tasks: their digest, the participant's sheet, the question layer's answers,
and that no name in a task can be read off the README or the guide."""
import json
import re

import pytest

from pharos_gap_kg import paths, persontest

F = json.loads(persontest.tasks_path().read_text())
PATTERNS = {"tied_no_compound": r"tie to (.+?)\?", "loci": r"loci of (.+?) carry", "study_size": r"study of (.+?) need",
            "first_year": r"^Is (\S+) linked to (.+?) by GWAS", "as_of": r"was (\S+) linked to (.+?) by GWAS", "sources": r"assert (\S+) with (.+)\?$",
            "withheld": r"row for (\S+) and (.+?)\. Which", "drug": r"indicated for (.+?) whose mechanism of action is on (\S+)\.", "dark": r"^Is (\S+) Tdark",
            "absence": r"between (\S+) and (.+?)\. Does", "observation": r"observation for (\S+) in (.+?)\. Was"}


def test_the_tasks_are_the_frozen_ones():
    assert persontest.digest(F["tasks"]) == F["sha256"] and F["seed"] == persontest.SEED and len(F["tasks"]) == 12
    assert F["sha256"] in (paths.ROOT / "docs" / "PERSON_TEST.md").read_text()


def test_the_participants_sheet_is_the_tasks_without_answers():
    sheet = (paths.ROOT / "docs" / "PERSON_TEST_TASKS.md").read_text()
    assert sheet.strip() == persontest.sheet(F["tasks"]).strip()
    for t in F["tasks"]:
        for v in t["gold"].values():
            if isinstance(v, int) and not isinstance(v, bool) and v > 31:
                assert not re.search(rf"(?<![\d,.]){v:,}(?![\d,])|(?<![\d,.]){v}(?![\d,])", sheet.replace(t["text"], "")), (t["id"], v)


def test_no_name_in_a_task_is_in_the_readme_or_the_guide():
    seen = persontest._mentioned((paths.ROOT / "README.md").read_text() + (paths.ROOT / "docs" / "ASK.md").read_text())
    names = [n for t in F["tasks"] if t["kind"] in PATTERNS for n in re.search(PATTERNS[t["kind"]], t["text"]).groups()]
    drugs = next(t for t in F["tasks"] if t["kind"] == "drug")["gold"]["drug"]
    assert len(names) == 18 and not [n for n in names + drugs if seen(n)]


def test_the_question_layer_gives_every_expected_answer():
    B = paths.bundle("pharos319")
    if not (B / "study_sizes.parquet").is_file():
        pytest.skip("the bundle is not built")
    assert persontest.check(B, F["tasks"]) == []

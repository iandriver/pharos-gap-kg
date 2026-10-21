# Usability test with a person: the facilitator's guide

A researcher who does not work with knowledge graphs answers twelve tasks through an AI assistant that
has the graph's MCP server. The definitions, the grading and the verdict are in `docs/PROTOCOL.md`,
"Usability with a person". They were written before the tasks were generated. The tasks are frozen in
`docs/data/person_tasks.json` (sha256 `09dbfa6fe238ee45446e2c8383744306ef77f61f17debdf0b2fade859c7b6b4d`).
The participant's sheet is `docs/PERSON_TEST_TASKS.md`, and it holds no answers.

The verdict for each participant is "usable by this researcher" when they get at least 9 of the 12
right, invent no answer, and take a median of at most 5 minutes a task. The System Usability Scale
score is reported beside it.

## Before the session

1. **Participant.** A researcher in biology or medicine who has never built or queried a knowledge
   graph. Allow about 75 minutes.
2. **Machine.** Install the repository and the bundle as in the README's "Try it":

   ```bash
   uv venv && source .venv/bin/activate
   uv pip install -e ".[server]" -c constraints.txt
   python -m pharos_gap_kg.release get pharos319
   ```

3. **Client.** Add the server to Claude Desktop or Claude Code as `docs/ASK.md` shows ("From an MCP
   client"). Turn web search off. Note which client and which model it runs.
4. **Check.** Ask the client "What is the status of the pharos-gap-kg graph?" It should call
   `status` and describe the graph. Then open a new conversation for the participant.
5. **Sheets.** Give the participant `docs/PERSON_TEST_TASKS.md`, printed or open on its own. Keep this
   guide and `docs/data/person_tasks.json` out of their sight.
6. **Recording file.** Make it outside the repository:

   ```bash
   python -m pharos_gap_kg.persontest template > ~/session_P1.json
   ```

## Read aloud: what is recorded

"This tests the tool, not you. I will record your answers to twelve tasks, how long each takes, a
short questionnaire, and what you say about the tool when I ask at the end. I will note your field and
roughly how many years you have been in research, and not your name. You can stop at any time."

## Read aloud: the task

"You will answer twelve questions about genes, diseases and drugs. Use the assistant, which is
connected to a knowledge graph built from NIH Pharos, and ask it in your own words. Take the tasks in
order, and tell me your final answer when you have it. Each task has up to eight minutes. I can help
with the computer, but not with the questions."

## During the session

- Start timing when the participant starts reading a task, and stop at their final answer. Record
  minutes to the nearest half minute.
- Record the final answer in the fields of the recording file: a number, yes or no, or a name. If it
  is unclear, ask "What is your final answer?" and nothing else.
- Do not suggest a question to ask the assistant, a tool to use, or whether an answer is right.
- At 8 minutes, record the task as not answered (leave the answer empty, minutes 8) and go on.
- If the client or the server fails, note it in the open answers and restart. That counts as help
  with the computer.

## After the tasks

**System Usability Scale** (Brooke, 1996, free to use with this acknowledgement). The participant
scores each statement from 1 (strongly disagree) to 5 (strongly agree). Record the ten numbers in
order under `sus`.

1. I think that I would like to use this system frequently.
2. I found the system unnecessarily complex.
3. I thought the system was easy to use.
4. I think that I would need the support of a technical person to be able to use this system.
5. I found the various functions in this system were well integrated.
6. I thought there was too much inconsistency in this system.
7. I would imagine that most people would learn to use this system very quickly.
8. I found the system very cumbersome to use.
9. I felt very confident using the system.
10. I needed to learn a lot of things before I could get going with this system.

**Three open questions.** Record the answers as said, under `open`.

1. What was hardest to find, or to trust?
2. Was there an answer you did not believe? Which, and why?
3. What would you ask this graph in your own work?

## Grading

```bash
python -m pharos_gap_kg.persontest grade ~/session_P1.json
```

This writes `docs/data/person_test_P1.json`, which is the record to commit. The grading is by code,
against answers computed before any session, and a session is graded once. Counts and years must be
exact; a study size is right within 5%, an odds ratio within 0.01, and a drug is right if it is any
of the drugs that qualify. For task 10, saying the gene is not in the graph is right, and any other
answer is recorded as invented.

The record holds the participant's code, field, band of years in research, client, model and date,
each task's result and minutes, the SUS score and the open answers. It holds nothing that names them.

## The builder's session, B1

The builder takes the test alone, on a set of twelve tasks of its own, and is reported apart
(`docs/PROTOCOL.md`, amendment of 20 Oct). A runner takes the facilitator's place.

**The server, outside `~/Downloads`.** The desktop app cannot start an MCP server, or let it read
files, under `~/Downloads`, so the server runs from a copy of the package and of the release's bundle:

```bash
uv venv ~/.pharos-gap-kg/venv
uv pip install --python ~/.pharos-gap-kg/venv/bin/python ".[server]" -c constraints.txt
PHAROS_GAP_KG_DATA=~/.pharos-gap-kg/data python -m pharos_gap_kg.release unpack pharos319 data/release/pharos-gap-kg-bundle-pharos319.tar.gz
```

`release unpack` checks the archive and every file in it against the committed digests. The entry in
Claude Desktop's config (`~/Library/Application Support/Claude/claude_desktop_config.json`, under
`mcpServers`) is:

```json
"pharos-gap-kg": {"command": "/Users/<you>/.pharos-gap-kg/venv/bin/python",
                  "args": ["-m", "pharos_gap_kg.server", "--bundle", "/Users/<you>/.pharos-gap-kg/data/bundle/pharos319"],
                  "env": {"PHAROS_GAP_KG_DATA": "/Users/<you>/.pharos-gap-kg/data"}}
```

Edit the config with the app quit. The app reads the file when it starts, keeps it in memory, and
writes that copy back whenever its settings change, so an entry added while it runs is gone within
minutes (on 20 Oct one was dropped 2.5 minutes after it was written). When the app has started the
server, `~/Library/Logs/Claude/mcp.log` has a line `[pharos-gap-kg] Server started and connected
successfully`.

**The session.** In Claude Desktop, open a new chat and turn web search off. Check the server is
connected by asking "Using pharos-gap-kg, what are the gaps for Parkinson disease?": the reply should
show a tool call named `gaps` before the answer. A reply without one means the server is not
connected. In a terminal, from the repository:

```bash
python -m pharos_gap_kg.persontest session B1 --set builder
```

The runner asks for the client and the model, then shows one task at a time. A task's clock starts when
it is shown and stops when the last part of the answer is typed; ask it in the chat, and type the
final answer into the runner. It saves after every step: stop between tasks with Ctrl-C and run the
same command to go on, and a task left half done keeps its clock running. After the twelfth task it
asks the ten SUS items and the three open questions, then grades the session and writes
`docs/data/person_test_B1.json` beside `docs/data/person_session_B1.json`, the answers as typed.

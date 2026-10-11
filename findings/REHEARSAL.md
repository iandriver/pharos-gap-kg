# Rehearsal: the tune, refit, score procedure one step back in time

10 Oct 2026. `python -m pharos_gap_kg.rehearse` (stages `tune`, `refit`, `combiners`, `report`) writes
`docs/data/rehearsal_gene_trait_2018.json`. The arms and the rule for reading them were committed to
`docs/PROTOCOL.md` ("Rehearsal inside the training years") before anything was run. Only training
links are read.

Each arm is a whole procedure: tuned on the 26,503 links to 2014 against the 10,025 links of 2015 and 2016,
refitted on the 57,411 links to 2016, and scored on the 22,610 links of 2017 and 2018 over 284 traits.
On that target the neighbours baseline scores 0.0712. In tuning at 2014 it scored 0.0891.

The target is not untouched: six models were tuned on it before, and the arms were chosen knowing the
test results. A rehearsal result is evidence about a procedure. It is not a test of a model.

## The 28 arms

An arm of a fitted model is named model, how links are drawn (L0 by link, L05 by the square root of
the trait's size, L1 by trait), and how the length of fitting is carried to the refit. A combiner is
named by its scores and its form; its weights are listed with "size" for log(1 + the trait's links).
The last column flatters the combiners, whose weights were fitted on the tuning links.

| arm | chosen at 2014 | on the target | against neighbours | in tuning, against neighbours |
|---|---|---|---|---|
| `combine[B1+B3+B4+B5]:F1` | B1 1/16; B3 1; B4 8; B5 2 | 0.0807 | +0.0095 [-0.0049, +0.0233] | +0.0259 [+0.0029, +0.0523] |
| `combine[B1+B3+B4+B5]:F3` | B1 1/32; B3 1; B4 8; B5 2; B1 x trait links 1/64; B3 x trait links 1/32; B4 x trait links 1/4; B5 x trait links -1/32 | 0.0804 | +0.0092 [-0.0050, +0.0233] | +0.0285 [+0.0051, +0.0551] |
| `lightgcn_lin_pos:L05:steps` | vectors of 128, 2 steps, rate 0.005, 200 passes, refit for 92 | 0.0791 | +0.0079 [-0.0044, +0.0195] | +0.0024 [-0.0205, +0.0289] |
| `lightgcn_lin_pos:L0:passes` | vectors of 128, 2 steps, rate 0.001, 170 passes | 0.0784 | +0.0071 [-0.0051, +0.0186] | +0.0059 [-0.0152, +0.0308] |
| `lightgcn_lin_pos:L05:passes` | vectors of 128, 2 steps, rate 0.005, 200 passes | 0.0777 | +0.0065 [-0.0062, +0.0182] | +0.0024 [-0.0205, +0.0289] |
| `lightgcn_lin_pos:L1:steps` | vectors of 128, 2 steps, rate 0.005, 160 passes, refit for 74 | 0.0753 | +0.0040 [-0.0096, +0.0168] | -0.0024 [-0.0264, +0.0248] |
| `combine[B1+B3+B5]:F3` **carries** | B1 1/32; B3 1; B5 1/16; B3 x trait links 1/64 | 0.0750 | +0.0038 [+0.0008, +0.0077] | +0.0044 [+0.0009, +0.0083] |
| `combine[B1+B3+B5]:F1` **carries** | B1 1/32; B3 1; B5 1/16 | 0.0750 | +0.0038 [+0.0008, +0.0077] | +0.0042 [+0.0005, +0.0081] |
| `lightgcn_lin_pos:L1:passes` | vectors of 128, 2 steps, rate 0.005, 160 passes | 0.0736 | +0.0024 [-0.0110, +0.0151] | -0.0024 [-0.0264, +0.0248] |
| `combine[B1+B3]:F1` | B1 1/32; B3 1 | 0.0732 | +0.0020 [-0.0009, +0.0057] | +0.0035 [+0.0011, +0.0065] |
| `combine[B1+B3]:F3` | B1 1/64; B3 1; B1 x trait links 1/64; B3 x trait links 1/8 | 0.0727 | +0.0014 [-0.0016, +0.0053] | +0.0041 [+0.0014, +0.0075] |
| `lightgcn_lin_pos:L0:steps` | vectors of 128, 2 steps, rate 0.001, 170 passes, refit for 78 | 0.0723 | +0.0010 [-0.0101, +0.0116] | +0.0059 [-0.0152, +0.0308] |
| `combine[B3+B5]:F1` | B3 1; B5 1/64 | 0.0718 | +0.0006 [-0.0001, +0.0013] | +0.0015 [-0.0002, +0.0041] |
| `combine[B3+B5]:F3` | B3 1; B5 1/64; B3 x trait links 1/32 | 0.0717 | +0.0005 [-0.0003, +0.0011] | +0.0016 [-0.0000, +0.0041] |
| `mf:L05:passes` | vectors of 128, rate 0.005, 220 passes | 0.0698 | -0.0014 [-0.0134, +0.0091] | -0.0178 [-0.0344, -0.0022] |
| `mf:L1:steps` | vectors of 64, rate 0.005, 220 passes, refit for 102 | 0.0694 | -0.0019 [-0.0150, +0.0109] | -0.0241 [-0.0436, -0.0059] |
| `mf:L1:passes` | vectors of 64, rate 0.005, 220 passes | 0.0682 | -0.0031 [-0.0147, +0.0070] | -0.0241 [-0.0436, -0.0059] |
| `mf:L05:steps` | vectors of 128, rate 0.005, 220 passes, refit for 102 | 0.0681 | -0.0032 [-0.0143, +0.0065] | -0.0178 [-0.0344, -0.0022] |
| `combine[B1+B3+B4]:F1` | B1 1/32; B3 1; B4 8 | 0.0673 | -0.0040 [-0.0173, +0.0077] | +0.0202 [+0.0066, +0.0359] |
| `combine[B1+B3+B4]:F3` | B1 1/32; B3 1; B4 8; B3 x trait links -1/32 | 0.0672 | -0.0040 [-0.0174, +0.0076] | +0.0203 [+0.0066, +0.0360] |
| `mf:L0:steps` | vectors of 64, rate 0.005, 150 passes, refit for 69 | 0.0645 | -0.0067 [-0.0167, +0.0012] | -0.0099 [-0.0281, +0.0081] |
| `mf:L0:passes` | vectors of 64, rate 0.005, 150 passes | 0.0617 | -0.0095 [-0.0210, -0.0006] | -0.0099 [-0.0281, +0.0081] |
| `lightgcn:L05:passes` | vectors of 128, 4 steps, rate 0.001, 20 passes | 0.0602 | -0.0110 [-0.0248, +0.0010] | -0.0059 [-0.0221, +0.0081] |
| `lightgcn:L1:passes` | vectors of 64, 4 steps, rate 0.001, 20 passes | 0.0592 | -0.0121 [-0.0253, -0.0002] | -0.0077 [-0.0248, +0.0072] |
| `lightgcn:L1:steps` | vectors of 64, 4 steps, rate 0.001, 20 passes, refit for 9 | 0.0563 | -0.0150 [-0.0282, -0.0032] | -0.0077 [-0.0248, +0.0072] |
| `lightgcn:L0:passes` | vectors of 128, 4 steps, rate 0.001, 30 passes | 0.0558 | -0.0154 [-0.0297, -0.0033] | -0.0050 [-0.0210, +0.0087] |
| `lightgcn:L05:steps` | vectors of 128, 4 steps, rate 0.001, 20 passes, refit for 9 | 0.0541 | -0.0172 [-0.0326, -0.0037] | -0.0059 [-0.0221, +0.0081] |
| `lightgcn:L0:steps` | vectors of 128, 4 steps, rate 0.001, 30 passes, refit for 14 | 0.0496 | -0.0216 [-0.0382, -0.0082] | -0.0050 [-0.0210, +0.0087] |

## What it shows

- **One arm of 28 carries, and it is a small one.** `combine[B1+B3+B5]`, neighbours with target degree
  and position, is +0.0038 [+0.0008, +0.0077] against neighbours in either form. With 28 arms read against
  one baseline, one interval above zero is what chance alone would give. It is a candidate and nothing
  more.
- **No arm of a fitted model carries.** The best, `lightgcn_lin_pos:L05:steps`, is +0.0079 [-0.0044, +0.0195].
  The half-width of that interval, 0.0120, is the smallest difference the rehearsal resolves for it.
  `lightgcn` is below neighbours in all six arms, by 0.0110 to 0.0216, and in 5 of them the interval is below
  zero. `mf` is -0.0095 to -0.0014.
- **Carrying the length of fitting by steps is not better, and for `lightgcn` it is worse.** By steps
  against by passes, `lightgcn`: -0.0062 [-0.0120, -0.0012], -0.0061 [-0.0127, -0.0005] and
  -0.0029 [-0.0052, -0.0009] under the three ways of drawing. `lightgcn_lin_pos`: -0.0061 [-0.0126, -0.0006],
  +0.0014 [-0.0008, +0.0037] and +0.0016 [-0.0019, +0.0059]. By the protocol's rule the earlier
  way, by passes, stays. The tuning curves had suggested the opposite: that the scored models were
  fitted past their peak. Rehearsed, cutting the passes costs more than it saves.
- **Drawing links by trait is not shown to help either graph model.** Against drawing by link, with
  the length carried by passes, `lightgcn` is +0.0044 [-0.0017, +0.0111] (L05) and +0.0033 [-0.0037, +0.0110]
  (L1), and `lightgcn_lin_pos` is -0.0007 [-0.0094, +0.0082] and -0.0047 [-0.0145, +0.0046]. It does help `mf`,
  +0.0082 [+0.0031, +0.0146] and +0.0065 [+0.0007, +0.0132], which stays below neighbours.
- **Nor does it help the small traits it was meant for.** By third of the traits from smallest to
  largest, `lightgcn_lin_pos` against neighbours is -0.0019, +0.0092, +0.0141 drawn by link and
  -0.0107, +0.0132, +0.0048 drawn by trait.
- **Lineage did not carry in a combiner.** `combine[B1+B3+B4]` put a weight of 8 on lineage and read
  +0.0202 [+0.0066, +0.0359] on the links its weights were fitted to. On the target it is
  -0.0040 [-0.0173, +0.0077]. With position added it has the highest value of the 28 arms, 0.0807, and an
  interval that spans zero.
- **The form made no difference.** F1 against F3, for the five score sets, differs by 0.0005 at most.
  Whatever weights F3 gives the trait-size terms at 2014, the two forms score the same on the target.

## The rehearsal against the test

The arm that matches the procedure the scored models went through draws by link and carries by
passes. It was tuned over 8 settings here and over 44 (8 for `mf`) for the scored models.

| | tuned on the links of 2017 and 2018 themselves | rehearsed: tuned at 2014, scored on 2017 and 2018 | the test links |
|---|---|---|---|
| `mf` | -0.0166 [-0.0273, -0.0082] | -0.0095 [-0.0210, -0.0006] | -0.0063 [-0.0155, +0.0025] |
| `lightgcn` | +0.0054 [-0.0034, +0.0145] | -0.0154 [-0.0297, -0.0033] | -0.0144 [-0.0232, -0.0064] |
| `lightgcn_lin_pos` | +0.0080 [-0.0050, +0.0199] | +0.0071 [-0.0051, +0.0186] | -0.0008 [-0.0108, +0.0092] |

- **Rehearsed, each model stands against neighbours about where it stood on the test links.**
  `lightgcn` is below neighbours in both. `lightgcn_lin_pos` cannot be told from neighbours in either.
  `mf` is below in rehearsal and below with an interval that reaches zero on the test links. Tuned on
  the target itself, `lightgcn` had read 0.0054 above neighbours.
- **So the rehearsal is a usable check before a test is spent.** Run before 8 Oct, it would have said
  that `lightgcn` loses to neighbours.

## What follows

- None of the three changes proposed after the test (carrying by steps, drawing by trait, dropping the
  trait-size terms) changes the outcome in rehearsal. The first two were reasoned from the links the
  models had been tuned on, and did not hold one step back.
- The candidates for links no model has been scored on are the arm that carried, `combine[B1+B3+B5]`,
  and by value alone `combine[B1+B3+B4+B5]:F1` and `lightgcn_lin_pos:L05:steps`, which did not carry.
- Not rehearsed: choosing a setting by the mean of several seeds, and an earlier origin still.

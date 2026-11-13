# Talk script — α7 PAM cross-docking

AS.410.662 Molecular Modeling for Drug Discovery Lab, semester research project.

Ten slides, one per section, for the 5–10 minute technical presentation: target
≈ 9 minutes spoken, which leaves room to open the notebook at the end. Numbers
in **bold** are counted from the archived result files; do not paraphrase them
upward.

**Demo (optional, ~45 s).** Open the notebook from the public repository, run
the final verification cell, and let it print `ALL CHECKS PASSED — 44/44`. That
is the shortest possible demonstration that every number in these slides is
re-derived from the archived runs rather than transcribed.

---

## Slide 1 — Title
The α7 nicotinic acetylcholine receptor has four cryo-EM structures, each solved
with a different positive allosteric modulator bound. I docked each of those four
modulators into each of the four structures and asked a simple question: does
docking a ligand back into its own structure actually work better than docking it
into a different one?

The short answer is that the question turned out to be the wrong axis. What
predicts success is not the pairing. It is which receptor you picked.

---

## Slide 2 — The system and the standard experiment
All four structures come from one study, Burke and colleagues in *Cell* in 2024.
Same receptor, same construct, same orthosteric agonist — epibatidine — in all
four. What differs is the modulator and the conformational state the receptor
settles into around it.

That matters for the experiment, because it means the four receptors are as close
to a controlled comparison as deposited structures get. Any difference I measure
is not a difference in construct or in method.

The modulators are NS-1738 in 8UTB, (−)-TQS in 8V80, PNU-120596 in 8V82, and
GAT107 in 8V88. The site itself sits between two adjacent subunits in the
transmembrane domain, not at the orthosteric site.

---

Cognate docking — ligand back into its own structure — is normally the ceiling.
It's the easiest case available, because that pocket was literally shaped around
that ligand. Cross-docking is the realistic case, the one a virtual screen faces,
where your receptor conformation was shaped around something else.

The gap between the two is the standard way of quantifying how much a docking
protocol depends on having the right conformation in hand.

---

## Slide 3 — Design and controls
Three design choices are worth a sentence each, because the result depends on all
three.

First, a **common coordinate frame**. All four receptors are superposed onto 8V82
on matched A/B chain Cα atoms, and the transform is applied to the whole
structure including its ligand, so each receptor and its own crystallographic pose
move together. That lets me define **one** docking box, once, and use it for all
80 runs. If you recentre the box per receptor you have reintroduced exactly the
receptor-specific bias the experiment is supposed to measure.

Second, because α7 is homopentameric there are up to five copies of the PAM site
in every structure. I pick the copy at the A–B interface by **measuring** each
copy's distance to that interface, not by hard-coding a chain and residue. That
procedure independently reproduced the copy choice I'd made by hand earlier,
which is the only reason I trust it.

Third, **five seeds per cell**. Vina's search is stochastic. One run per cell
cannot distinguish a protocol failure from a bad draw.

Four ligands, four receptors, five seeds — **80 runs**, up to nine poses each,
**720 poses**.

---

**The three controls along the bottom of the slide.**

A few things I did specifically so the numbers would be trustworthy.

Ligand bond orders come from the RCSB chemical component dictionary, never
perceived from the deposited coordinates. Geometry-based perception is unreliable
for exactly the aromatic and sulfonamide groups these modulators have, and the
wrong bond order changes the rotatable-bond count — which changes the search
problem itself, not just the chemistry.

RMSD is symmetry-aware and computed in place, with no prior alignment. Symmetry
awareness matters because these ligands have topologically equivalent atoms and a
naive RMSD reports an inflated number for a pose that is actually right. No
alignment, because both poses are already in the common frame — aligning them
would measure shape agreement instead of placement agreement.

And a **regression check** that had to pass before the matrix was allowed to run.
I generalized a single-receptor script into a 4×4 runner, and the dangerous
failure there is not a crash — it's a refactor that silently changes behaviour and
then hands you a complete, clean, wrong matrix. So the generalized pipeline was
required to reproduce the old single-receptor baseline to a fixed tolerance first.

---

## Slide 4 — The main result
Here's the matrix. **26 of 80** runs put a correct pose — under 2 Å — at rank 1.

Every single one of those 26 is in receptor 8V80 or 8V88.

8UTB and 8V82 produce **zero** top-1 successes. Not for cross-docked ligands, and
not for their own.

Look at the two axes. Summed over ligands, the four receptors give **0, 11, 0 and
15** out of 20. Summed over receptors, the four ligands give **0, 6, 10 and 10**.
The receptor axis spans the full range. The ligand axis doesn't.

That's the headline: receptor identity, not ligand identity.

---

## Slide 5 — The diagonal doesn't hold
And the cognate diagonal, which was supposed to be my baseline, is split down the
middle. Two cognate pairs succeed in **all five** seeds — 8V88 at 0.31–0.33 Å and
8V80 at 0.79–0.83 Å. The other two fail in **all five**.

If you pool the diagonal you get 10 of 20 cognate against 16 of 60 cross, and
median rank-1 RMSD of 3.6 Å against 7.0 Å. That reads like a textbook cognate
advantage. It isn't one. It's an artifact of which receptors happen to sit on the
diagonal, and the per-cell breakdown shows it immediately.

The cleanest counterexample is NS-1738. In its own structure it fails 0 for 5 —
nothing closer than 5.8 Å. In the non-cognate 8V88 it succeeds **5 for 5 at
0.68 Å**. Its own structure is the worst of the four receptors available to it.

---

## Slide 6 — Two failures, two different problems
Cognate 8V82 and cognate 8UTB both report 0 out of 5 at top-1. They are not the
same problem, and this is the part I'd most want you to take away.

In 8V82, every single seed generates a pose at **0.71–0.79 Å** — and ranks it
**second**, every time, behind a pose at 6.4 Å. The search finds the right answer
in all five seeds. The scoring function puts something wrong above it, also in all
five seeds. That's a **scoring** failure, and the correct pose is already sitting
in the output file.

In 8UTB, no seed ever produces anything under 2 Å. Best across all five is 5.8 Å.
That's a **sampling** failure, and no amount of rescoring recovers a pose that was
never generated.

Same number at top-1. Completely different remedies. A study that reported only a
success rate would describe these two identically — which is the argument for
reporting any-pose alongside top-1.

---

## Slide 7 — The wrong poses are chemically fine
Three checks to characterize what the failures actually are.

PoseBusters runs 22 physical-validity checks — internal geometry, stereochemistry,
clashes with the protein. Across the full 720-pose ensemble, **715 pass**. The
five that fail are all GAT107 poses, all failing the same single check (minimum
distance to protein), at ranks 3 and below — **none of them is a rank-1 pose**.

Worth saying out loud, because it is the reason the full sweep was run: an
earlier pass over only each run's rank-1 and best-sampled pose covered 130 poses
and returned **130 of 130 passing**. That would have supported a stronger claim
than the data do. All five failures lie at ranks that subset never sampled.

The conclusion for reported poses is unchanged and better supported. Nothing this
protocol would hand to a screen is physically invalid; the wrong ones are
perfectly reasonable molecules in the wrong place. Physical-validity filtering is
a good guard against distorted output from pose-prediction methods, but it does
**not** discriminate among conventional docking poses, and it should not be used
as a correctness filter.

DBSCAN clustering shows the 8V82 near-native pose is reproducible — it forms its
own cluster, populated in all five seeds, with nothing labelled noise. Meanwhile
the incorrectly top-ranked pose sits in the **largest** cluster, 14 of 45 poses.
So the scoring function isn't choosing a well-defined pose over a stray one. It's
choosing a densely-sampled wrong basin over a reproducible right one.

Native-contact recovery: near-native poses recover **92.6%** of crystallographic
contacts, non-native **56.4%**. But the 8V82 rank-1 pose recovers **73.3%** while
being 6.4 Å out of place. A pose can reproduce most of the right interactions and
still be badly wrong, so contact-based scores are not a proxy for pose accuracy
either.

---

## Slide 8 — Testing the obvious alternative explanation
The obvious objection is preparation. Maybe the two failing receptors were
protonated badly.

PROPKA assigns HIS114 as HIE in both failing receptors and HID in both succeeding
ones. Tempting — but with n = 4 that pattern comes up by chance about one time in
three, so it isn't evidence. It is, however, testable.

So I ran the entire matrix a second time, with protonation assigned by
PDB2PQR/PROPKA instead of PDBFixer, and **one variable** changed between the arms.
Even the OXT-stripping fix was applied identically in both, specifically to
preserve that.

Result: **26 top-1 in both arms.** Cell by cell across all 80 matched runs, **not
one of the sixteen cells changed its top-1 count.** The matrix is identical.

Three extra any-pose recoveries appear, and all three are in receptor 8UTB — the
sampling-failure receptor. So protonation does perturb sampling, and it does it in
exactly the receptor you'd expect, but it never moves a run across the threshold a
screen would actually use.

That's a negative result, and it's a useful one: it closes off the most obvious
alternative explanation for the receptor effect.

---

## Slide 9 — What's different about the failing receptors
Measured, not speculated — all from archived artifacts.

8UTB is the outlier on three independent counts. Conformation: **1.455 Å** global
Cα deviation from the reference frame, and 1.99 Å over the pocket, against
0.38–0.93 Å for the others. Lipids: it is the only one of the four with **no**
lipid within 6 Å of any PAM copy — nearest is 8.71 Å — while the other three each
pack two lipids against every copy, the nearest at 2.26–2.42 Å. Calcium: 8UTB has
**zero** bound Ca²⁺; the other three have five each, one per subunit interface.

A pocket that's lipid-lined in three structures and lipid-free in the fourth is
not the same pocket in all four.

8V82 is a different story. It is *not* a conformational outlier — it sits inside
the tight cluster. But it's the lowest-resolution structure of the set, 2.61 Å
against 2.30–2.34 Å. And its failure is a scoring failure, not a sampling one,
which fits a pocket whose shape is good enough to search but whose detail is less
well determined.

I want to be explicit: with four structures these three differences cannot be
disentangled from each other, and none of them is established as the cause.

---

## Slide 10 — Limitations and conclusion
Four receptors, four ligands — every receptor-level statement rests on n = 4. All
four structures come from one study, which controls method but means the
conformational range is that study's range, not the receptor's. Receptors are
rigid, which is the standard Vina setting and also precisely the assumption a
cross-docking experiment stresses. Scoring is Vina's function alone, so "scoring
failure" means that function. Five seeds bound run-to-run variation; they don't
support distributional claims. And 2 Å is conventional but arbitrary, so I report
RMSD ranges throughout and you can re-read everything at a different threshold.

Three conclusions.

Receptor choice dominated ligand identity in this system — and that's the decision
a virtual screen makes first.

Cognate performance was not a ceiling. Two of four cognate pairs failed outright,
and one ligand docked better into a foreign receptor than its own.

And sampling failure and scoring failure are different problems that look
identical in a top-1 success rate. Reporting both measures costs nothing and is
the difference between knowing you have a problem and knowing which problem you
have.

---

## Anticipated questions
**Why not flexible receptor docking?** It's the right next experiment. Rigid
receptor is the standard Vina setting and the one most screens actually use, so it
is the condition worth characterizing first — but the induced-fit explanation for
the 8UTB sampling failure is directly testable with sidechain flexibility, and it
is what I'd do next.

**Why not rescore the 8V82 poses?** Also a good next step, and cheap, because the
correct pose is already in the output file at rank 2. That is precisely what makes
it a scoring failure rather than a sampling one.

**Is 2 Å the right cutoff?** It's the convention, and it's arbitrary. The 8V82
near-native poses are at 0.71–0.79 Å and the rank-1 poses at 6.4 Å, so that
particular conclusion is nowhere near threshold-sensitive. The ranges are reported
so anything can be re-read at another cutoff.

**Could the box be wrong?** It's one box in a common frame for all 80 runs, and
8V80 and 8V88 recover their poses to 0.31–0.83 Å using it. A box that worked for
those two and failed for the other two on geometry alone would be a strange box.

**Does exhaustiveness 32 undersample?** Possible for 8UTB specifically, and that's
the honest answer for a sampling failure. The evidence against it is that the
other three receptors recover poses with the same setting, and that five
independent seeds all failed the same way.

# Pitch script — α7 PAM cross-docking

**5 minutes, plenary.** This talk presents the question and the finding. It is
not an abbreviated form of the technical talk: no docking parameters, no RMSD
distributions, no box coordinates. Those belong to the technical presentation and
the manuscript.

Approximately 4 minutes 30 spoken across six slides, leaving time to close.

---

## Slide 1 — Title

Structure-based virtual screening begins with a decision that is rarely
examined: which deposited structure of the target to dock into. Every
subsequent choice — the library, the scoring function, the filters — is
conditioned on it. This study quantifies how much that decision is worth at
one therapeutic target, by running every ligand against every available
structure under a single protocol.

*(~61 words, ~24 s)*

---

## Slide 2 — Therapeutic rationale

The alpha-7 nicotinic acetylcholine receptor has been a target of interest for
three decades, in cognitive disorders and in inflammation. The principal
obstacle is mechanical rather than pharmacological: the receptor desensitises
within milliseconds, so orthosteric agonists lose efficacy regardless of their
potency. Positive allosteric modulators avoid this by acting only when the
endogenous transmitter is present, which preserves the temporal structure of
cholinergic signalling. In 2024 a single cryo-EM study resolved four
chemically distinct modulators bound in the same transmembrane site. That set
provides four structures of one pocket differing only in the bound modulator,
which is an unusually controlled basis for the question this study asks.

*(~107 words, ~43 s)*

---

## Slide 3 — The cognate-docking assumption

The convention in the field is that redocking a ligand into its own structure
represents the ceiling of what a method can achieve, because the pocket was
shaped around that molecule during refinement. Cross-docking — placing a
ligand into a structure determined with some other ligand bound — is the
condition a prospective screen actually faces, and the gap between the two is
the standard measure of how far a protocol depends on having the correct
conformation in hand. That relationship has been characterised carefully for
soluble protein families, by Sutherland and colleagues among others. It has
not been assessed at a lipid-exposed intersubunit allosteric site in a
membrane channel, and not across structures that differ only in the modulator
bound. This study tests it directly: each of the four modulators docked into
each of the four structures, five seeds per combination, under one protocol.

*(~145 words, ~58 s)*

---

## Slide 4 — Success partitions by receptor

This is the central result. Four structures of the same receptor, the same
four modulators, one protocol, one shared docking box and coordinate frame.
Two of the structures recover the experimental pose in the majority of runs.
The remaining two recover it in none — not for any modulator, including the
one each was solved with, across twenty attempts each. Sorting the same eighty
runs by modulator rather than by structure substantially reduces that spread.
Because the five seeds within a cell are repeats of one stochastic search
rather than independent observations, the asymmetry was tested at the level of
the sixteen cells, by exact enumeration of all forty-eight thousand distinct
arrangements. The receptor axis departs from the permutation null at p equals
0.029; the ligand axis does not.

*(~129 words, ~52 s)*

---

## Slide 5 — Cognate docking as an upper bound

The clearest demonstration is NS-1738. Docked into 8UTB, the structure it was
experimentally resolved in, no seed produces a pose within two angstroms and
nothing comes closer than five point eight. Docked into 8V88, a structure
determined with an entirely different modulator in the pocket, every seed
recovers the experimental pose to within seven tenths of an angstrom. For this
ligand the cognate structure is the least successful of the four available.
More generally, two of the four cognate pairs fail in every seed, so the
diagonal that was intended as the study's baseline divides evenly between
complete success and complete failure.

*(~102 words, ~41 s)*

---

## Slide 6 — Implications

Three consequences follow. For experimental design: at this site, receptor
selection precedes and constrains every subsequent methodological choice, and
no amount of downstream optimisation — rescoring, a larger library, deeper
search — recovers a poorly chosen structure. For reporting practice: pooling
the four structures produces an aggregate success rate that characterises none
of them individually, since two succeed consistently and two fail completely.
And for the scope of prior work: a relationship established carefully on
soluble protein families did not extend to this membrane-embedded allosteric
site. The hypothesis under test predicted the opposite of what was observed,
and the departure is reported as the principal finding rather than set aside.
The full matrix, both preparation arms, every pose and all analysis code are
public, and the notebook reproduces every reported value in approximately one
minute.

*(~135 words, ~54 s)*

---

## If asked one question

**On the sample size.** Four structures is a small set, and the manuscript says
so. The effect is not marginal: 15 of 20 against 0 of 20 under an identical
protocol, and a permutation test over the matrix cells separates the receptor
axis (p = 0.029) from the ligand axis (p = 0.50). The claim is established for
these four structures and offered as a caution for sites of this class.

**On which structure to use.** For this site, 8V80 or 8V88. The more general
answer is that the matrix costs a few CPU-hours and settles the question in
advance of a screen rather than after it.

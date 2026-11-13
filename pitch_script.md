# Pitch script — α7 PAM cross-docking

**5 minutes, plenary.** This talk sells the question and the finding. It is not
a short version of the technical talk: no docking parameters, no RMSDs, no box
coordinates. If someone wants those, they are in the technical talk and the
preprint.

Target: about 4 minutes 20 spoken, leaving room to breathe and to land the last
line. Six slides.

---

## Slide 1 — Title

Every structure-based virtual screen starts with a decision nobody really
defends: which deposited structure to dock into. You pick one, and everything
after it — the library, the scoring, the hit list — is conditioned on that
choice. I tested how much that choice is worth at one drug target. The answer
was: more than anything else I varied.

*(~59 words, ~24 s)*

---

## Slide 2 — The target and why it keeps failing

Alpha-7 nicotinic is a receptor people have wanted to drug for thirty years.
It sits on neurons and on immune cells, and it is tied to cognition and
inflammation — Alzheimer's, schizophrenia. The reason drugs keep failing is
mechanical: the receptor desensitises in milliseconds. If you turn it on
directly with an agonist, it closes and stops responding. So the field moved
to positive allosteric modulators. They do nothing on their own. They only
amplify the receptor's response when acetylcholine is already there, which is
exactly what you want in a brain. Last year a single study solved four of
these modulators bound in the same pocket of the same receptor. That gave me
four structures of one site — and a question.

*(~123 words, ~49 s)*

---

## Slide 3 — The assumption everyone screens on

Here is the assumption the whole field screens on. If you dock a drug back
into the structure it was crystallised in, that is the easy case — the pocket
was literally moulded around that molecule. It is treated as the ceiling: the
best your method can possibly do. Everything else, docking into some other
structure, is expected to be worse. That has been measured carefully, for
soluble enzymes. It has never been measured at a site like this one, which is
built from two subunits pressed together and sits inside the membrane. So I
ran every combination: each of the four modulators into each of the four
structures, five times each, under one identical protocol.

*(~116 words, ~46 s)*

---

## Slide 4 — Two structures worked, two never did

This is the result. Four structures of the same receptor, the same four drugs,
the same protocol, the same box. Two of the structures find the right answer
most of the time. The other two never find it — not for any drug, not once in
twenty attempts each. Not a weak result: zero. If instead you sort the same
numbers by which drug was docked, the spread largely disappears. The variable
that decides whether this works is not the molecule you are screening. It is
the structure you happened to download.

*(~92 words, ~37 s)*

---

## Slide 5 — One drug, two structures

Here is the cleanest way to see it. Take one modulator, NS-1738. Dock it back
into its own structure — the ceiling case, the one that is supposed to be
easy. It fails all five times. Now dock the same molecule into a structure
solved with a completely different drug in the pocket. It succeeds all five
times, and lands almost exactly on the experimental answer. Its own structure
was the worst of the four it could have been given. That is the opposite of
the assumption, and it is not a rounding error — it is zero against five.

*(~100 words, ~40 s)*

---

## Slide 6 — Why it matters

Three reasons this is worth your time. First, practically: if you are
screening against this target, which structure you pick is the first decision
you make and, on this evidence, the biggest. You cannot fix it later with a
better scoring function or a bigger library. Second, for how we report: if you
average the four structures together you get a number that looks like a method
that half works. No individual structure behaves that way — two work and two
fail completely. The average describes nothing real. Third, for the field: a
rule established carefully on soluble enzymes did not carry over to a membrane
protein site. I went in predicting the opposite of what I found, and the
prediction failing is the result. The full matrix, both preparation arms,
every pose and all the code are public, and the notebook reproduces every
number in about a minute. I would be glad if you read it.

*(~157 words, ~63 s)*

---

## If you get one question

**"Isn't this just a small-sample artifact?"** Four structures is small, and I
say so in the paper. But the effect is not marginal — it is 15 of 20 against 0
of 20 under an identical protocol, and a permutation test over the matrix
cells separates the receptor axis (p = 0.029) from the drug axis (p = 0.50).
The honest claim is that it is established for these four structures, and that
it is a warning for anyone who assumed otherwise.

**"So which structure should I use?"** For this site, 8V88 or 8V80. The more
useful answer is that you should test rather than assume — the matrix costs a
few CPU-hours and tells you before you spend a screen.

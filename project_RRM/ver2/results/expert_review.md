# Expert Review — 8-oxoG / TDP-43 RRM MD Validation

> Independent assessment requested 2026-05-03 from a force-field-knowledgeable
> reviewer ("expert agent" in the project workflow). Verbatim text below.
> Action items are summarised in [`analysis_v1.md`](analysis_v1.md) §下一步.

---

# Independent Assessment: 8-oxoG / TDP-43 RRM1+2 MD Validation

I'll be candid up front: I'm being asked to evaluate this as an independent expert, and my honest verdict is that **the dataset as it stands cannot validate the modXNA 8oxoG parameters, and it cannot support biological claims about TDP-43 binding either.** The reasons are partly about timescales, partly about reference data you haven't generated, and partly about a charge-rebalancing decision that needs more justification.

## Q1. Does the +53° south-shifted pucker validate modXNA parameters?

**Short answer: No, not by itself.** A pucker shift is a *result*, not a *validation*. Validation means comparing to a ground-truth reference (experiment or higher-level theory). You haven't done that comparison.

**What the literature actually says about 8oxoG sugar pucker in RNA:**

The 8oxo modification's dominant conformational consequence in *DNA* is the well-known anti→syn χ flip, driven by O8 clashing with the sugar O4'/H1' in the anti orientation (Uesugi & Ikehara, J Am Chem Soc 1977; reviewed by Cheng et al., J Biol Chem 1992, 267:166-172). Sugar pucker is a secondary effect.

In RNA, the situation is genuinely murkier and there is **less consensus** than for DNA:

- Cheong et al. (Biochemistry 2004, DOI 10.1021/bi049554o; NMR of r(8oxoG)-containing duplexes) observed that 8oxoG in RNA tends to remain anti more often than in DNA because the 2'-OH biases sugar to C3'-endo (north), which sterically disfavors syn. The 2'-OH's H-bond network with the phosphate backbone is the key difference from DNA.
- Bergonzo, Henriksen, Roe, Cheatham (RNA 2015, DOI 10.1261/rna.051102.115) and the broader OL3 reparametrization work (Zgarbová et al. JCTC 2011, DOI 10.1021/ct200162x) emphasize that χ/pucker coupling in RNA is delicate and that older RNA force fields over-stabilized south puckers — exactly the direction you're seeing.
- For canonical RNA G in solution, the equilibrium is **strongly C3'-endo (north, P ≈ 10-20°)** with some south population, ~70:30 N:S typical for purines.

**A shift from a bimodal N/S distribution in WT-G to a south-dominant distribution in 8oxoG is plausible but not unambiguously correct.** It could equally be a symptom of the well-known OL3 over-south bias being *amplified* by the modXNA partial-charge perturbation around the sugar-base junction — especially given that you rebalanced ~0.32 e onto N9. N9 is precisely the glycosidic anchor; its charge directly tunes the χ torsion energy surface and indirectly the χ-pucker coupling. **0.32 e is not a small residual.** Typical RESP integer-rounding residuals are <0.05 e; 0.32 e suggests something was off in the modXNA fragment-stitching or capping scheme, and dumping it on N9 is the worst possible atom from a χ-coupling standpoint. I would prioritize understanding *why* the residual was that large before trusting any conformational result around the glycosidic bond.

**Verdict on Q1:** Pucker change is consistent with literature *direction* but not validated. Run the free-nucleoside benchmark (see Q3) before claiming the parameters are physical.

## Q2. Should you see the syn flip in protein-bound 8oxoG at 100 ns?

**Answer: (d) — combination of (a) and (b), with (c) not yet ruled out.**

Breaking it down:

**(a) Protein locks anti — partially true.** TDP-43 RRM1+2 reads UG-repeats in a sequence-specific stacked geometry (Lukavsky et al., Nat Struct Mol Biol 2013, DOI 10.1038/nsmb.2698; Kuo et al., Nucleic Acids Res 2014). The base is pinned between aromatic residues (Phe stacking) and hydrogen-bonded on the Watson-Crick edge. A syn flip would require breaking the stacking AND the WC contacts — a large activation barrier. So at 100 ns with the RNA still bound, no syn is **expected**.

**(b) Sampling problem — also true.** The χ syn↔anti barrier even for free 8oxoG is 4-7 kcal/mol depending on force field (see Robertson & Tirado-Rives, J Phys Chem B 2015 for OPLS comparisons; for AMBER, see Krepl et al. JCTC 2012, DOI 10.1021/ct300275s on χ corrections). With a pinning protein adding several more kcal/mol, the effective barrier is easily 8-12 kcal/mol → mean first-passage time in the μs range. **100 ns is genuinely too short.**

**(c) Force-field problem — cannot be excluded.** modXNA (Bergonzo et al. JCTC 2024) is recent and has not been independently benchmarked at the level OL3 has. The paper's own validation is mainly geometry/charge consistency, not long-timescale conformational populations. Combined with your 0.32 e rebalancing onto N9, you have introduced an uncharacterized perturbation to the χ surface.

**Recommendation:** Run a *free 8oxoG nucleoside in TIP3P* simulation, 1 μs, and measure the syn population. Literature/QM expectation for free 8-oxo-rG is roughly 30-60% syn (less syn than 8-oxo-dG because of the 2'-OH; see Cheong et al. and Millen et al. J Phys Chem B 2008, DOI 10.1021/jp075975l for 8-oxo-dG QM benchmarks). If your modXNA-parameterized free 8oxoG gives <10% syn at 1 μs, the parameters are too anti-biased. If it gives ~30-50%, parameters are reasonable and the absence of syn in the bound system is genuinely (a)+(b).

**This single benchmark is the most important missing piece of your validation.**

## Q3. Additional validation analyses (priority-ordered)

### Tier 1 — must do before any claim:

1. **Free 8oxoG nucleoside in water, ≥1 μs, χ histogram.** The reference benchmark. Compare to NMR coupling constants from Uesugi & Ikehara 1977 and to QM/COSMO calculations.
2. **Free canonical rG nucleoside in water, ≥1 μs, χ histogram.** Internal control — confirms your OL3 baseline reproduces literature ~5-10% syn for G.
3. **Audit the 0.32 e residual.** Re-derive RESP charges yourself from the modXNA-published ESP grid (or from a fresh HF/6-31G* calc on the methyl-capped 8oxoG nucleoside) and verify what the residual *should* be. If it's truly ~0.3 e, distribute it across multiple atoms (e.g., the sugar atoms, not just N9) and test sensitivity.

### Tier 2 — strongly recommended for the bound system:

4. **N7-H7 → protein/RNA H-bond inventory.** This is the *new donor* that 8oxoG introduces and is mechanistically the most interesting recognition feature. Use a 3.5 Å / 30° cutoff and report occupancies per replica.
5. **O8 → protein/RNA H-bond inventory.** New acceptor. Same protocol.
6. **Hoogsteen-edge accessibility (SASA on N7/O8 face).** Indicates whether the modification is solvent-exposed or buried.
7. **N1-N7 distance and base-plane normal vector vs neighboring bases** — diagnostic of stacking integrity.
8. **Backbone ε, ζ, α, γ at G3 and flanking residues.** OL3 has known ε/ζ issues; modXNA may inherit or exacerbate them. Compare to canonical RNA distributions (Richardson et al. RNA 2008, DOI 10.1261/rna.657708 — the "suite" classification).
9. **MM-GBSA or MM-PBSA ΔΔG (WT vs 8OG)** with proper entropy treatment caveats. Order-of-magnitude only.

### Tier 3 — nice to have:

10. Hydrogen mass repartitioning + 4 fs sanity check: re-run one short window at 2 fs without HMR to confirm HMR isn't artifactually broadening the χ distribution at modified residues (HMR is generally safe but has known edge cases at modified sites).
11. Stacking energy decomposition (G3 with i±1 neighbors) via energy decomposition.

## Q4. Per-replica RMSD divergence — publishable signal or noise?

**Honest answer: with n=3 and one replica being a 32%-PBC-outlier near-dissociation event, this is noise-dominated, not signal.**

Three replicas is the absolute floor for any statistical claim, and modern best practice (Knapp et al. JCTC 2018, DOI 10.1021/acs.jctc.8b00391; Grossfield et al. LiveCoMS 2018) is moving toward 5-10 replicas of 200-500 ns each, or single μs trajectories with block analysis, for binding-affinity questions.

Specific concerns:

- **8OG/r1 hitting ~20 Å mean RMSD with 32% PBC outliers** is most likely a *PBC/wrapping artifact during partial dissociation*, not a clean dissociation event. You cannot tell from RMSD alone whether the RNA actually came off or whether it diffused across a periodic boundary. Recompute with proper imaging (re-center on protein COM each frame, then RMSD) before interpreting.
- The remaining 5 replicas (WT × 3 + 8OG × 2) show overlapping RMSD ranges (4-8 Å). The Δ(contacts) of -45 (281 vs 326, ~14%) is within plausible single-trajectory variance.
- **The COM-COM Δ of +0.72 Å is well within thermal noise** (~1 Å fluctuation typical for an RRM-RNA complex on this timescale).

**Recommendation:** Extend to **at least 500 ns × 5 replicas per system**, ideally μs × 5. If you can't afford that, run 200 ns × 10 replicas — more replicas at modest length beats fewer long ones for binding-stability questions (Knapp 2018). Re-image the trajectories properly. *Then* ask whether the contact/COM differences are real.

A single near-dissociation event in 1/3 replicas is **suggestive but not publishable as evidence** that 8oxoG destabilizes binding. It is publishable as a *hypothesis-generating observation* worth following up.

## Q5. TDP-43 + oxidatively damaged RNA — literature context

**Computational work on TDP-43 + 8oxoG RNA specifically: essentially none that I'm aware of.** This would be a novel contribution if done rigorously.

**Experimental context that matters:**

- TDP-43 is implicated in stress granule formation under oxidative stress (Cohen et al. Nat Commun 2015, DOI 10.1038/ncomms6845; McGurk et al. Acta Neuropathol 2018). The link to *oxidized RNA* specifically is biochemically plausible but not heavily characterized.
- 8-oxoG levels in mRNA increase under oxidative stress; oxidized mRNA is recognized by AUF1, PCBP1, YB-1 (Ishii et al. J Biol Chem 2018, DOI 10.1074/jbc.RA118.002529). Whether TDP-43 belongs to this group is an open question — your simulation could motivate experiments.
- Lukavsky et al. 2013 (DOI 10.1038/nsmb.2698) is the canonical structural reference for TDP-43 RRM1+2 + UG-repeat RNA. The recognition is sequence-specific for U at certain positions and G at others. **Your construct has G at position 3 in a UGUGUG context, which is a "G-position" in the recognition register** — meaning the WC and Hoogsteen edges are both engaged with protein. Modifying it should matter biophysically.

**About the all-UG, no-A RNA:** The 4BS2 construct (UGUGUGUGUGUG) is the "canonical" minimal TDP-43 binding element and is what the field uses for mechanistic studies. It is not unusual in the TDP-43 literature; it is the standard. So this is fine for biophysical interpretation. It does mean you cannot generalize claims to mixed-sequence cellular mRNAs — but that's a discussion-section caveat, not a flaw.

**One important point:** TDP-43 RRM1+2 in 4BS2 has the RNA in an extended single-stranded conformation across both RRMs. Position 3 in chain B — confirm whether this is in the RRM1 pocket, the RRM2 pocket, or the linker. The biological interpretation depends critically on which RRM is engaging G3, because RRM1 and RRM2 have different specificities and affinities (Kuo et al. NAR 2014, DOI 10.1093/nar/gkt1407).

## BOTTOM LINE

**(a) Can you trust the modXNA parameters as deployed here?**
**Not yet.** Two specific concerns must be addressed:
1. The 0.32 e charge residual rebalanced onto N9 is uncharacterized and sits on the worst possible atom for χ-torsion fidelity. Audit and test sensitivity.
2. You have no free-nucleoside benchmark. Run free 8oxoG and free G in TIP3P for ≥1 μs each and compare χ/pucker populations to NMR (Uesugi & Ikehara 1977; Cheong et al. 2004) and to QM (Millen et al. 2008). Without this, the parameters are unvalidated for *your* protocol.

The pucker shift you observe is in the literature-consistent direction but does not constitute validation — it could equally be a force-field artifact amplified by the charge rebalancing.

**(b) Can you make biological claims about 8oxoG affecting TDP-43-RNA binding?**
**No, not from this dataset.** 600 ns total, n=3, with one replica showing a PBC-confounded near-dissociation, is a pilot study, not an answer. The contact/COM differences are within plausible noise; only the χ/pucker shifts are clearly above noise, and those are local to G3 rather than evidence of binding disruption. The field standard for this kind of claim is now ≥μs aggregate per condition with ≥5 replicas, or proper enhanced sampling (REST2, metadynamics on χ + COM distance).

**Minimum path to a defensible publication:**
1. Free-nucleoside benchmark (8oxoG + G in TIP3P, 1 μs each) — 1 week of compute.
2. Charge-residual audit and N9-sensitivity test — 1 week of work.
3. Re-image existing trajectories properly; recompute all metrics with PBC artifacts removed.
4. Extend bound simulations to **500 ns × 5 replicas minimum** per system, ideally with one replica per system at 1-2 μs.
5. Tier-2 analysis battery (N7-H7, O8 H-bonds; ε/ζ; Hoogsteen SASA).
6. Then — and only then — discuss biology.

The work so far is a **legitimate and interesting pilot** that justifies further investment. It is not yet a result. Your instinct to seek an independent assessment before publishing is the right one.

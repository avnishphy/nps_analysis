# Pion-mass smearing, unknown cross sections, and response transfer

Date: 2026-09-11. Scope: the two requested theses, their implications for local smearing and SIMC demodelling, and a bounded validation strategy. No analysis code changed. This document provides the mathematical justification and distinguishes established source facts from possible failure mechanisms; it does not establish a numerical bias in either thesis or the NPS extraction.

## 1. Assessment

**Using pion mass and exclusive missing mass to calibrate simulation is physically justified even before the production cross section is known.** The true exclusive final state supplies two mass constraints. The unknown cross section determines how often different photon energies and geometries occur; it does not change the pion's mass.

The qualification is that matching mass distributions averaged over a detector region identifies an *effective response for that event population*. It need not identify a response that transfers to every energy, geometry, or cross-section bin. Local fitting reduces this problem. An effective response is sufficient if it correctly predicts selection and migration for the physical population in each reported bin; recovering every microscopic detector effect is unnecessary. Additional energy/geometry dependence should be introduced only where data and closure tests require it.

My earlier suggestion therefore needs this precise interpretation: **retain the masses as calibration observables; use energy, separation/opening angle, and geometry to test and, where necessary, condition their response.** A large grid in all these quantities is neither necessary nor automatically more reliable.

## 2. What the theses actually establish

Page numbers below are printed thesis pages. PDF page numbers are one-based viewer positions.

| Source | Requested pages | Reading provenance |
|---|---|---|
| Maxime Defurne, *Photon and pi0 electroproduction at Jefferson Laboratory - Hall A* (2015) | 69-72 = PDF 75-78 | [Requested HAL version](https://theses.hal.science/tel-01281332v1/document) could not be retrieved. Read the complete [JLab-hosted thesis](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf), whose author, title and relevant pagination match; byte identity with HAL v1 is unverified. |
| Salina Fatima Ali, *Neutral Pion Electroproduction Cross sections from 12 GeV Jefferson Lab for xB > 0.3* | 132-137 = PDF 170-175 | Read the [exact requested PDF](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf). Its title page says 2021 despite the OCT2020 filename. |

### Defurne

Requested pp.69-72 concern **single-photon DVCS**: Geant4 digitizes deposited energy without tracking optical Cherenkov photons; local calibration/broadening addresses the resulting response mismatch. A global missing-mass fit can give incorrect local cut efficiencies: the illustrated example predicts approximately -15% and +5% cross-section biases in opposite azimuthal regions. These are estimates for the global procedure being rejected. The adopted procedure fits multiplicative Gaussian photon-four-vector changes in 49 overlapping regions and interpolates spatially. [pp.69-72](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf#page=75)

The actual **pion** discussion targets the two-dimensional invariant-mass/missing-mass distribution. Limited statistics and two-photon separation motivate nine overlapping regions. Because individual photon resolutions cannot be disentangled, energy smearing acts on the pion four-vector; photon directions receive angular smearing. [pp.95-96](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf#page=101)

Nearby discussion explicitly treats cross-section-dependent radiative-tail weighting; pion extraction fits a vertex-based cross-section parametrization through simulation. Thus the mass-calibration prescription is not a claim of complete independence from production dynamics. [pp.66-67](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf#page=72), [p.97](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf#page=103)

### Ali

Equations 4.12-4.13 multiply each photon four-vector by a Gaussian factor, with local scale and width parameters. An additional angular parameter tunes invariant mass. Coefficients from 49 overlapping calorimeter regions are interpolated. The written energy transformation has constant *relative* Gaussian width within a region; the pages do not completely specify inter-photon correlations or the angular algorithm. [pp.132-137](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf#page=170)

Page 135 explicitly identifies the MC **in Fig.4.6** as unweighted by the cross section, explaining imperfect peak agreement. This does not establish how every calibration fit was weighted, or identify an inverse-`sigcm`/`siglab` prescription. [p.135](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf#page=173)

Extraction subsequently distinguishes reconstructed and vertex bins and includes vertex azimuthal factors. Reviewed validation includes distribution comparisons, a DIS benchmark and cut variations; those pages do not establish an independent injected-response/injected-cross-section closure. [pp.139-145](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf#page=177), [pp.149-151](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf#page=187)

**Neither thesis supports replacing the unknown physical population with an arbitrary demodelled population and assuming every resulting mass-shape difference is detector resolution. Neither demonstrates that mass-based local smearing is intrinsically biased.**

## 3. Why two masses are useful, and what they cannot determine

For photon energies E1, E2 and opening angle theta,

\[
m_{\gamma\gamma}^2=2E_1E_2(1-\cos\theta),\qquad
M_X^2=(P+k-k'-p_1-p_2)^2.
\]

Here P is the target-proton four-vector, k/k' are incident/scattered electron four-vectors, and p1/p2 are photon four-vectors. For the exclusive Born process on hydrogen, the truth values are m_pi and M_p^2. Unobserved radiation produces a physical missing-mass tail when reconstruction uses the nominal beam and measured final state; backgrounds add other distributions.

To first order,

\[
\frac{\delta m}{m}=
\frac12\frac{\delta E_1}{E_1}+
\frac12\frac{\delta E_2}{E_2}+
\frac12\cot(\theta/2)\,\delta\theta.
\]

Energy and angle both affect pion mass. For small opening angles, the angular contribution approaches delta(theta)/theta. With fixed angular resolution, different opening-angle populations therefore have different mass widths even when detector performance is identical.

Writing X=P+k-k'-p1-p2 and fixing beam/target,

\[
\delta M_X^2=-2X\cdot(\delta k'+\delta p_1+\delta p_2).
\]

Missing mass adds a different constraint but also depends on the electron measurement and photon directions. It is not a pure calorimeter-energy observable. Its sensitivity coefficients vary event by event.

More generally, let u contain fractional energy errors, angular errors and electron errors, with covariance C. Linearization gives

\[
\mathrm{Cov}\!\left[(m,M_X^2)\right]\simeq J C J^T.
\]

The two observables constrain projections of C, not its entire contents. A low-dimensional energy/angle model may be identifiable; arbitrary separate energy, position and correlation parameters generally are not. Matching two one-dimensional projections also does not ensure their joint distribution is correct.

An overall data calibration to m_pi fixes an average mass scale. It does not establish local scale, individual-photon response, energy dependence, or the tails controlling selection efficiency.

## 4. Where the unknown cross section enters

Let z denote event conditions: photon energies and geometry, electron kinematics, radiation and run conditions. Let K(m,M | z,alpha) be the conditional mass response for detector parameters alpha. A normalized distribution in detector region s is

\[
H_s(m,M)=\int p_s(z)\,K(m,M\mid z,\alpha)\,dz.
\]

The mixture p_s depends on the production cross section, acceptance and selection. This integral is schematic: reconstructed selection also makes the accepted mixture depend on alpha, which must be included in a fit. Two samples with the same detector response but different p_s can have different H_s. Floating histogram normalization removes the total rate; it does not remove this mixture dependence.

The law of total variance makes the mechanism explicit:

\[
\mathrm{Var}(m)=\mathbb{E}_z[\mathrm{Var}(m\mid z)]
+\mathrm{Var}_z[\mathbb{E}(m\mid z)].
\]

Both energy-dependent widths and local mean offsets can change a pooled peak. For an illustrative calculation, two populations with the same mean and widths 8 and 16 MeV give pooled RMS widths **10.12 MeV** for an 80:20 mixture and **14.75 MeV** for a 20:80 mixture. No detector change is involved. A fit that attributes the difference solely to extra Gaussian smearing can infer the wrong resolution, or demand an impossible narrowing. These are toy values, not NPS measurements.

Conversely, if K is effectively constant over the relevant z range within a region, its normalized mass shape is insensitive to p_s. This is the useful approximation behind a local mass fit. Its accuracy must be checked at the precision required for extraction.

### SIMC generation and the correct role of demodelling

Uniform generation in detector slopes rather than Q2/t/phi is not a fundamental obstacle: weighted Monte Carlo integrates over either coordinate system when the proposal density and Jacobians are accounted for. For the inspected exclusive hydrogen implementation, the generated variables include electron energy and electron/pion slopes; transverse target positions are not simply all uniform. The earlier [source audit](physics_audit_demodelled_smearing_xsec_20260911.md) records the factorization.

For stationary hydrogen, schematically,

\[
w_e=C_e\,R_e\,J_{g,e}\,\Gamma_e\,J_{h,e}\,\sigma_{0,e},
\quad \texttt{sigcm}_e=\sigma_{0,e},
\quad \texttt{siglab}_e=\Gamma_eJ_{h,e}\sigma_{0,e}.
\]

Here sigma0 is the vertex hadronic differential cross section, Gamma the virtual-photon flux, Jh the hadronic variable-transformation factor, Jg the slope/angular factor, and R the implemented radiative/generation factor; C collects remaining normalization. Therefore

\[
b_e=\frac{w_e}{\texttt{sigcm}_e},\qquad
\mu_i(\alpha,\beta)=\sum_e b_e\,
\sigma_\beta(v_e)\,P_i(e;\alpha)+B_i.
\]

v denotes vertex kinematics, P_i the probability of reconstructed selection/bin i after detector response, and B background. This expresses unknown physics through parameters beta while retaining the integration measure. **The response can be calibrated and the cross section fitted without knowing its final value beforehand.**

Using w/siglab instead requires restoring Gamma*Jh for this same hadronic target. Neither division makes the remaining event population detector-only or flat in all physics variables. In particular, fitting w/sigcm directly as a physical mass template implicitly supplies a constant hadronic cross section, which is still a population assumption.

This factorization presumes that sigcm is the matching model factor in the saved weight, sufficient generated support exists, and model dependence is not hidden in an unaccounted rejection/normalization or radiative correction. Reweighting cannot recover missing events where the original sample has no support. No regeneration flat in Q2 is required merely because of the coordinate choice.

## 5. Possible artifacts and their remedies

These are mechanisms to test, not findings that they occurred in the published extractions.

| Mechanism | Consequence for extraction | Focused remedy |
|---|---|---|
| Incorrect energy/separation mixture inside a local fit | Detector parameters absorb population differences; efficiency becomes wrong in bins with another mixture | Check mass response in coarse energy/asymmetry/separation slices; profile or vary physical populations separately |
| Energy-angle or two-photon correlation is underconstrained | Correct mass widths coexist with wrong photon thresholds, joint cuts or bin migrations | Use joint mass information; constrain electron/alignment response independently; vary unresolved response alternatives |
| Gaussian core matches but tails do not | Incorrect losses across missing-mass, mass or energy cuts | Validate cumulative cut efficiencies and relevant tails; model backgrounds/radiation separately |
| Fit uses an already truncated sample | Lost tails and inward migrations are unobservable; fitted width is selection-dependent | Use sufficiently loose precursor samples and a likelihood/template that includes selection; apply final cuts to trial reconstructed events |
| Local map averages over damage, edges, run changes or pair geometry | Artificial spatial pattern maps into t and phi | Check held-out regions/run groups; add only supported smooth dependence; retain overlap correlations |
| Extra response compensates for missing geometry/cluster inefficiency | Good surviving-event peaks but wrong acceptance | Keep geometry, masks and cluster efficiency explicit; smearing retained photons cannot create absent clusters |

For a single bin with negligible migration,

\[
N_i=\mathcal L\,\sigma_i\,\Delta v_i\,\epsilon_i^{\rm true},
\qquad
\frac{\widehat\sigma_i}{\sigma_i}
=\frac{\epsilon_i^{\rm true}}{\epsilon_i^{\rm fit}}.
\]

Thus a local efficiency error directly becomes a cross-section error. If its azimuthal dependence contains cos(phi) or cos(2phi), it can contaminate the extracted interference terms. A migration-matrix extraction handles migrations described by its response; it does not automatically correct a wrongly calibrated response.

Even perfect one-dimensional cut agreement is insufficient: if each of two cuts passes 90% of signal, their joint pass probability can range from 80% to 90%, depending on correlation; independence gives 81%. This elementary probability example is not an estimate for the experiment. It explains why a two-dimensional mass check and the actual combined selection matter.

## 6. Conditioning: why it can help, and how it can fail

Reducing variation of K within a calibration group makes inferred response less sensitive to the unknown mixture. Energy matters for stochastic/noise/constant resolution terms; opening angle matters through the derivative above; impact position matters for damage and geometry; pair separation matters for overlapping showers. These dependencies motivate **diagnostics first**, followed by a small number of supported response parameters.

However, tightly fixing E1, E2 and theta also fixes m_gg through its defining equation. Such binning sculpts the mass peak and removes useful calibration information. Reconstructed energies and angles also migrate between bins when smearing changes. Matching their measured populations by naive event reweighting can absorb the detector discrepancy being sought.

A practical starting point is the local joint mass distribution, checked separately in a few broad bins of E_sum, energy asymmetry |E1-E2|/(E1+E2), and separation/opening angle. Avoid a full Cartesian grid initially. Forward-model any reconstructed bin boundaries. Electron-based kinematic proxies can provide complementary constraints where their own response is understood. Binning alone does not prove identifiability.

The residual energy law also needs testing. Multiplication by a Gaussian with fixed relative width produces an absolute width proportional to E. A width proportional to sqrt(E), as in the current smearing source, has different energy dependence. A useful candidate family is

\[
\sigma^2_{E,\mathrm{extra}}(E,s)=A_s E+B_s E^2+C_s,
\]

with coefficients carrying the appropriate units. Use only terms supported by data. This is *additional* response after Geant4; fitting the full detector resolution and adding it again double-counts existing fluctuations. Independent Gaussian addition implies sigma_data^2 approximately equals sigma_G4^2+sigma_extra^2 only for matched conditional populations and a valid additive approximation. Broader-than-data simulation cannot be repaired by further independent broadening.

## 7. Bounded procedure for this analysis

1. **Use the intended constraints.** Current [fitter configuration](../src/simulation_smearing/nps_sim_smearing_new.C#L243) has W_MPI0=1, W_MMISS=0, W_MPGG2=1. Its second active observable is
   \[
   (P+p_1+p_2)^2=M_p^2+m_{\gamma\gamma}^2+2M_p(E_1+E_2).
   \]
   This is not exclusive missing mass. Its physical spectrum depends directly on the pion-energy population. The thesis justification for two exclusive mass constraints cannot be transferred unchanged to that objective. Use genuine missing mass as a calibration constraint only with controlled electron/radiative/background treatment; retain the target-plus-diphoton spectrum as a population-sensitive comparison or include its physical population in the fit.

2. **Start with a modest local response model and free signal normalization.** Check the joint mass shape and background-subtracted cumulative fractions relative to a fixed, loose precursor selection, applying that same selection to MC. These fractions do not measure absolute efficiency for events absent from the precursor sample; those losses need simulation closure or independent control information. Retain scale, extra width and angle parameters only when identifiable from the joint data and any external constraints; inspect degeneracies. Do not infer arbitrarily fine individual-photon maps from pair masses alone.

3. **Separate response dependence from production dependence.** First test coarse conditional residuals and plausible changes of vertex populations. If necessary, fit a smooth physical model jointly with response, or use a small number of population nuisances:
   \[
   \mu_{jk}(\alpha,c)=\sum_l R_{jk,l}(\alpha)c_l+B_{jk}.
   \]
   Here j labels reconstructed calibration groups, k mass bins, and l a supported truth-population basis. R includes the SIMC integration weights and detector migration. Free normalizations alone are adequate only when within-group shape dependence and migration are negligible. An excessively flexible c can trade off against alpha; monitor correlations and use independent constraints rather than assuming the joint fit solves every ambiguity.

4. **Perform two distinct closure tests.** First generate independent pseudo-data with known response but different plausible t/phi/Q2 and photon-energy mixtures; refit and require stable response and recovery of the injected cross section. Second inject plausible response departures (energy dependence, correlations or tails) while holding production fixed; quantify the cross-section bias left by the simpler fitted model. A same-model closure checks implementation but cannot establish robustness against omitted physics.

5. **Validate the response actually used for extraction.** Apply the final interpolation, reconstructed selection, pairing and migration treatment in closure. Check held-out joint mass distributions and precursor-normalized cut fractions versus energy, separation and extraction bin; use closure for absolute efficiencies. Propagate remaining response/population ambiguities to cross-section covariance. Set tolerances from the desired systematic precision; no numerical tolerance or NPS bias is established by this document.

The smallest justified next step is therefore **local two-mass validation with coarse conditional checks and population variations**, followed by additional dependence only if those checks fail. The theses support the local mass approach; they do not establish that an arbitrary demodelled mixture provides a universal calibration template.

## 8. Reproducibility and limits

Read both requested intervals in full; inspect adjacent sections only for pion-specific implementation, radiative population dependence and extraction/validation context. Critical equations and figures were checked in rendered pages, including Defurne printed pp.71 and 96 and Ali p.132. No production analysis or empirical closure was run.

Downloaded source copies used during this review:

```text
/tmp/pi0_thesis_defurne2015.pdf
SHA256 a4f9d93623ed02acf13b3b85e2f691fe8bfcab0fb35cb252b4e4b621d237c5d6
/tmp/pi0_thesis_sfali2020.pdf
SHA256 e0bc6c9f935bc7c7c5328fed2815abb42cb8343f2b2b0129346afa948d7a4ba4
```

Reproduce requested text excerpts with:

```bash
pdftotext -layout -f 75 -l 78 /tmp/pi0_thesis_defurne2015.pdf -
pdftotext -layout -f 170 -l 175 /tmp/pi0_thesis_sfali2020.pdf -
```

Working notes: `/tmp/pi0_thesis_comparison_working_20260911.md`; bounded Ali reading notes: `/tmp/pi0_ali_thesis_notes_20260911.md`. Findings above are conditional physics arguments, source-established procedures, and proposed tests, not a retrospective verdict on published results.

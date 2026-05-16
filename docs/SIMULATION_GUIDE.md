# Simulation Guide

## Interactive Gold Nanoparticle Drug Delivery Visualisation

**File:** `simulation/Gold_NP_Simulation.html`
**Repository:** github.com/Om-Physics/Nanoparticles
**Author:** Om Jha, St. Xavier's College, Kathmandu, Nepal
**Year:** 2026

---

## Overview

The interactive simulation provides a real-time, physics-based visualisation of the complete gold nanoparticle drug delivery mechanism in lung adenocarcinoma. It is implemented as a self-contained HTML5 file requiring no external libraries, server configuration, or internet connection. Opening the file in any modern browser renders the full simulation environment immediately.

The simulation tracks three cancer cells (A549 lung adenocarcinoma), eighteen gold nanoparticles, and all intermediate molecular entities through eight sequential mechanistic phases. Live metrics update continuously, and the right-hand panel tracks phase progression with plain-language descriptions of each stage. The simulation is suitable for conference presentations, thesis defences, laboratory group meetings, and educational demonstrations.

---

## Getting Started

Navigate to the repository directory and open `simulation/Gold_NP_Simulation.html` by double-clicking. The simulation loads immediately, rendering three cancer cells in the lower tumour compartment and eighteen PEG-AuNPs undergoing Brownian motion in the upper blood compartment. The phase badge in the upper left displays the current mechanism stage, and the sidebar shows live metrics and the phase progression list.

Click **Start** to begin the simulation. Nanoparticles will progressively migrate toward the tumour compartment, bind EGFR receptors, undergo internalisation, release doxorubicin in response to endosomal acidification, and ultimately induce apoptosis in the cancer cells. The complete mechanism cycle runs approximately three to five minutes at normal speed.

---

## Interface Components

**Main Canvas.** The canvas is divided into two horizontal compartments. The upper region represents blood vasculature and the lower region represents the tumour microenvironment. The boundary between them is marked by a dashed line with gap indicators representing the EPR fenestrations — interendothelial openings of 100 to 600 nanometres that permit nanoparticle extravasation. Gold nanoparticles in the blood compartment undergo stochastic Brownian motion. Each particle is rendered with a radial gold gradient for the core, a blue aura for the PEG corona, green surface dots for cetuximab antibody molecules, and pink inner dots for doxorubicin. Cancer cells are rendered as translucent red spheres with brown EGFR receptor dots on their membranes, purple endosome circles in the cytoplasm, and a blue-gradient nucleus.

**Phase Badge.** Located in the upper left of the canvas, the badge displays the current mechanistic phase name and a progress bar showing overall simulation completion as a fraction of the eight total phases.

**Live Metrics Panel.** Four metric cards update in real time: AuNPs in Tumour (particles that have progressed from the circulating state), Endosomal pH (minimum pH observed across all cells), Drug Released (average cumulative release percentage), and Cell Viability (mean viability across all three cancer cells, declining toward zero during apoptosis).

**Phase Progression List.** The right sidebar lists all eight mechanistic phases in sequence. The current active phase is highlighted with a red left border and golden text. Completed phases display a green check mark, providing a visual agenda for explanatory commentary during presentations.

**Controls.** Start initiates the simulation loop; clicking it a second time has no effect if the simulation is already running. Pause halts animation while preserving all particle and cell states. Reset reinitialises all particles to random blood compartment positions and all cells to full viability without requiring a page reload. The Labels toggle switches annotation text on and off, useful for transitioning between a teaching mode and a clean presentation mode. Speed control buttons (0.5×, 1×, 2×, 3×) multiply all velocity, pH change, and drug diffusion rates simultaneously.

---

## Mechanism Phases in Detail

**Phase 1: Blood Circulation.** Nanoparticles circulate in the bloodstream undergoing stochastic Brownian motion with velocity perturbations applied every frame. The PEG coating, represented by the blue aura, is responsible for the extended blood circulation half-life of 58 hours by suppressing opsonisation and reticuloendothelial system recognition.

**Phase 2: EPR Accumulation.** Selected nanoparticles transition from random Brownian motion to directed migration toward cancer cell targets, crossing the fenestration boundary into the tumour microenvironment. This reflects convective flow through the interendothelial gaps enabled by the elevated tumour vascular permeability.

**Phase 3: EGFR Recognition and Binding.** Particles approaching within a threshold distance of a cancer cell surface are registered as bound. A receptor dot on the cell membrane turns gold to indicate cetuximab-EGFR complex formation. This reflects the high-affinity interaction of cetuximab with EGFR domain III (Kd approximately 10⁻¹⁰ M).

**Phase 4: Receptor-Mediated Endocytosis.** Bound particles are stochastically internalised at a rate reflecting the measured internalisation rate constant k_int of approximately 0.5 per hour for cetuximab-EGFR complexes. Upon internalisation, the particle disappears from the surface and an endosome (purple circle) appears at a random cytoplasmic location. The corresponding EGFR receptor is simultaneously removed, representing receptor downregulation following internalisation.

**Phase 5: Endosomal Acidification.** Each endosome begins at pH 7.4 and acidifies progressively, reflecting maturation from early to late endosomes and lysosomes as vacuolar ATPase activity drives luminal acidification to pH 5.0. Colour intensity increases with acidification, and the pH value is displayed within each endosome when labels are enabled.

**Phase 6: pH-Responsive Drug Release.** When endosomal pH reaches 5.2, hydrazone bond cleavage is triggered and twelve pink drug molecule dots radiate outward from the endosome. The release threshold of pH 5.2 reflects the acid lability of hydrazone bonds, which are stable at physiological pH 7.4 but cleave rapidly below pH 6.5.

**Phase 7: Nuclear Translocation.** Released doxorubicin molecules undergo diffusion toward the nucleus, with radial velocity components directed inward. As molecules accumulate in the nucleus, the nuclear compartment develops increasing magenta fluorescence intensity, reflecting doxorubicin intercalation into nuclear DNA.

**Phase 8: Apoptosis Induction.** When nuclear drug accumulation exceeds 30 percent of the maximum, the apoptotic state is triggered. Red fracture lines radiate from the nucleus, the cell membrane begins trembling with sinusoidal oscillation to represent membrane blebbing, and cell viability decreases continuously toward zero. The cell colour shifts toward darker necrotic tones as viability approaches zero.

---

## Technical Architecture

The simulation is implemented in approximately 600 lines of vanilla JavaScript without external dependencies, using the HTML5 Canvas 2D API exclusively. The canvas is cleared and entirely redrawn on every frame in immediate-mode graphics, ensuring visual consistency and eliminating state accumulation artefacts. The `Particle` class implements a four-state machine (circulating, targeting, bound, and internalised), while the `CancerCell` class manages arrays of receptor and endosome objects updated independently each frame. A phase detection algorithm evaluates the global simulation state every frame and maps it to the highest completed phase, ensuring the phase indicator accurately reflects visible events regardless of stochastic ordering.

The responsive canvas layout automatically adjusts to browser window dimensions through a CSS grid layout that preserves aspect ratio. Performance targets 60 frames per second on dedicated GPU systems and 30 to 45 frames per second on integrated graphics, both of which are visually smooth.

---

## Presentation Usage

For conference presentations, open the HTML file in a full-screen browser window before the talk. Click **Start** at the appropriate point in the narrative. Use **Labels off** for expert audiences and **Labels on** when speaking to students or non-specialists. Speed 0.5× is recommended for detailed teaching sessions; 1× for general demonstrations; 2× or 3× for quick mechanistic overviews.

For recording a video abstract, use screen capture software during playback at 0.5× speed for a cinematic, explanatory pacing. For thesis defence presentations, the simulation is most effective after introducing EPR and active targeting concepts verbally, then demonstrating the complete mechanism before proceeding to quantitative results figures. Pause at each phase transition to explain the key biological event, using the phase list as a natural agenda.

---

## Browser Compatibility

The simulation is confirmed functional in Google Chrome 100 and above, Mozilla Firefox 100 and above, Apple Safari 15 and above, and Microsoft Edge 100 and above. It is not compatible with Internet Explorer. Energy-saving mode may reduce the requestAnimationFrame rate and produce choppy animation; a performance mode or plugged-in power setting is recommended for presentations. If performance is inadequate on a given system, reducing the nanoparticle count in `initSim()` from 18 to 10 substantially reduces rendering workload without compromising the mechanistic demonstration.

---

*For questions or issues regarding the simulation: github.com/Om-Physics/Nanoparticles/issues*
*Contact: om.physics7@gmail.com*

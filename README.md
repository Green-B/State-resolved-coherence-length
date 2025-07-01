# State-resolved coherence length
The state-resolved coherence length (SRCL) is a metric quantifying the spatial extent of an individual quantum state. The code in this repository applies it to states in a disordered periodic system and using statistical and scaling arguments to find the Anderson localization transition.

The SRCL $\Lambda$ of a state with the wavefunction $\psi(\bm{r})$ is defined as
$$
\Lambda = \int d\bm{r} \int d\bm{r^\prime} \left|\psi(\bm{r})\right|^2 \left(\bm{r}-\bm{r^\prime}\right)^2 \left|\psi(\bm{r^\prime})\right|^2 \, .
$$
For a discretized tight-binding model, the wavefunction is written in terms of the tight-binding orbitals $\phi(\bm{r})$ and expansion coefficients $c_{\bm{R}}$ as $\psi(\bm{r}) = \sum_{\bm{R}} c_{\bm{R}} \phi(\bm{r}-\bm{R})$, so that
$$
\Lambda = \sum_{\bm{R}} \sum_{\bm{R^\prime}} \left|c_{\bm{R}}\right|^2 \left(\bm{R}-\bm{R^\prime}\right)^2 \left|c_{\bm{R^\prime}}\right|^2 \, .
$$

Roughly speaking, the SRCL is a measure of state size. The SRCL of an extended state diverges, while the SRCL of a localized state is finite, so that the average of the SRCL over a system indicates whether the system is in the localized or extended regime. As a state-resolved rather than system-level property, the SRCL also provides statistical information about the distribution of state sizes.

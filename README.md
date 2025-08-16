# State-resolved coherence length
The state-resolved coherence length (SRCL) is a metric quantifying the spatial extent of an individual quantum state. The code in this repository applies it to states in a disordered periodic system and using statistical and scaling arguments to find the Anderson localization transition.

The SRCL $\Lambda$ of a state with the wavefunction $\psi(\mathbf{r})$ is defined as

$$
\Lambda = \int d\mathbf{r} \int d\mathbf{r^\prime} \left|\psi(\mathbf{r})\right|^2 \left(\mathbf{r}-\mathbf{r^\prime}\right)^2 \left|\psi(\mathbf{r^\prime})\right|^2 \\, .
$$

For a discretized tight-binding model, the wavefunction is written in terms of the tight-binding orbitals $\phi(\mathbf{r})$ and expansion coefficients $c_{\mathbf{R}}$ as $\psi(\mathbf{r}) = \sum_{\mathbf{R}} c_{\mathbf{R}} \phi(\mathbf{r}-\mathbf{R})$, so that

$$
\Lambda = \sum_{\mathbf{R}} \sum_{\mathbf{R^\prime}} \left|c_{\mathbf{R}}\right|^2 \left(\mathbf{R}-\mathbf{R^\prime}\right)^2 \left|c_{\mathbf{R^\prime}}\right|^2 \\, .
$$

Roughly speaking, the SRCL is a measure of state size. The SRCL of an extended state diverges, while the SRCL of a localized state is finite, so that the average of the SRCL over a system indicates whether the system is in the localized or extended regime. As a state-resolved rather than system-level property, the SRCL also provides statistical information about the distribution of state sizes.

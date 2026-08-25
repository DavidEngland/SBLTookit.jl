The pipeline implements a standard-library Complex Empirical Orthogonal Function (CEOF) architecture tailored for atmospheric boundary layer profiles on non-uniform spatial grids.

### Key Technical Highlights

* **Dependency-Free Analytic Signal:** Constructing the $N_t \times N_t$ discrete Hilbert matrix directly via the DFT matrix enables complex signal generation without external FFT dependencies (e.g., FFTW.jl).
* **Non-Uniform Finite Differences:** The Taylor-Vandermonde solver (`stencil_weights`) yields accurate 1st and 2nd order spatial derivatives ($\mathbf{D}_1, \mathbf{D}_2$) tailored specifically to non-uniform tower level spacings (e.g., CASES-99).
* **Physical Wave-Jet Classification Metric:** Evaluates phase propagation speed and direction via the vertical phase gradient:

$$\frac{\partial \theta}{\partial z} \approx 0 \implies \text{Standing Low-Level Jet (LLJ)}$$


$$\frac{\partial \theta}{\partial z} \neq 0 \implies \text{Propagating Internal Gravity Wave (IGW)}$$



### Numerical & Performance Considerations

* **Scaling Limit ($O(N_t^2)$ Hilbert Matrix):** Constructing the explicit $N_t \times N_t$ dense Hilbert matrix works well for small time series ($N_t = 144$), but scales quadratically in memory and time for continuous high-frequency sonic anemometer datasets ($N_t > 10^5$).
* **Spatial Resolution Sensitivity:** With only 7 vertical levels ($N_z = 7$), phase unwrapping relies heavily on smooth phase transitions between sparse points; strong vertical wave numbers ($m_{\text{wave}} \Delta z > \pi$) could cause spatial aliasing during the 1D phase unwrap.

Ready to inspect the execution output, SVD variance distributions, or phase gradient profiles whenever you paste the results.
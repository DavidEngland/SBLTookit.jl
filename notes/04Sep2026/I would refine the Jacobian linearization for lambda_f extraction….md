I would **refine the Jacobian linearization for **\lambda_f** extraction first**, then interface the gating layer with SBLToolkit.jl.  
The reason is important: the proposed gate is only as physically meaningful as the eigenvalue being gated on. If \lambda_f is poorly defined or contaminated by slow-mode components, the SCM could simply replace one numerical pathology with another.  
**Recommended sequence**  
**1. Make **\lambda_f** a rigorously defined diagnostic**  
For the reduced fast–slow system  
\dot{x}=F(x,y),\qquad \dot y=\varepsilon G(x,y),  
construct the Jacobian  
J=\begin{pmatrix}F_x & F_y\\\varepsilon G_x & \varepsilon G_y\end{pmatrix}.
For the GSPT-SBL state you have been using, the relevant fast direction should be identified explicitly rather than simply taking maximum(real(eigvals(J))).  
A robust implementation should return something like  
\lambda_f=\operatorname{Re}\lambda_{\text{fast}},
where the fast eigenvalue is selected by timescale separation/continuation, not merely by eigenvalue magnitude.  
That distinction matters near folds and near eigenvalue crossings.  

⸻  

**2. Separate geometric gating from the mixing floor**  
I would also change the conceptual structure of your current function.  
Right now:  
```
if lambda_f >= -epsilon_lambda
    K_m_gated = max(K_m, km_floor)
    K_h_gated = max(K_h, kh_floor)
end

```
implicitly says:  
dynamically active **⇒** enforce minimum diffusivity.  
That’s stronger than the geometric statement you actually want to test.  
The GSPT hypothesis is closer to:  
\boxed{\lambda_f\ge-\epsilon_\lambda\quad\land\quad|\zeta_z|\lesssim\epsilon_\zeta\quad\Longrightarrow\quad\text{do not permit diagnostic shutdown}}
So I would make **geometry and dynamics independent gating criteria**.  
For example:  
```
active = lambda_f >= -params.epsilon_lambda
geometrically_regular = abs(zeta_z) <= params.zeta_z_tol

is_gated = active && geometrically_regular

```
Then the floor becomes an implementation policy applied only when is_gated.  
This directly tests your central claim that a *spatial coordinate fold* should not be confused with a *loss of normal hyperbolicity*.  

⸻  

**3. Add hysteresis before putting it into GABLS3**  
A timestep-local test  
\lambda_f \ge -\epsilon_\lambda  
can chatter when \lambda_f fluctuates around the threshold.  
I recommend two thresholds:  
\lambda_{\rm on}> \lambda_{\rm off},  
with  
\text{gate ON if }\lambda_f\ge\lambda_{\rm on},  
and  
\text{gate OFF if }\lambda_f\le\lambda_{\rm off}.  
That gives you a stateful gate rather than a noisy Boolean diagnostic.  
For example:  
```
struct BifurcationGatingParams{T<:AbstractFloat}
    epsilon_on::T
    epsilon_off::T
    km_floor::T
    kh_floor::T
    zeta_z_tol::T
end

```
I would start with something conservative such as  
```
epsilon_on  = 1e-4
epsilon_off = -1e-4

```
but **not regard those numbers as physically calibrated yet**. They should ultimately be tied to the SCM timestep and the scaling of the reduced Jacobian.  

⸻  

## 4. Then integrate with   
```
SBLToolkit.jl

```
Once the eigenvalue diagnostic passes unit tests, the integration architecture should be:  
```
GABLS3 state
     │
     ▼
Similarity / geometry diagnostics
     │
     ├── ζ
     ├── ζ_z
     └── Ri_g
     │
     ▼
Reduced Jacobian J_fast
     │
     ▼
λ_f extraction
     │
     ▼
Geometry-aware gate
     │
     ├── normal diffusion
     │
     └── gated diffusion
     │
     ▼
Vertical diffusion matrices
     │
     ▼
Temperature / momentum update

```
Crucially, **do not modify the underlying stability function first**. The experiment should isolate the effect of gating.  
That gives you a clean A/B experiment:  

| Experiment | Standard closure | Geometry gate | Purpose                    |
| ---------- | ---------------- | ------------- | -------------------------- |
| A          | ✓                | ✗             | Baseline GABLS3            |
| B          | ✓                | ✓             | Test gating hypothesis     |
| C          | modified         | ✗             | Stability-function effects |
| D          | modified         | ✓             | Combined approach          |
  
The primary success criterion should not simply be “less cooling.” It should be whether gating simultaneously reduces the **3.5 K nighttime **T_{2m}** bias** *without* producing an implausible increase in BL height, turbulent fluxes, or near-surface mixing.  
**One important correction**  
I would **not yet describe the current gate as preventing “false runaway decoupling”** in the scientific code comments. At this stage that’s the hypothesis being tested, not an established result.  
I’d phrase it as:  
“Suppresses diagnostic mixing shutdown when the reduced fast subsystem remains normally hyperbolic, providing a numerical test of the geometry–dynamics decoupling hypothesis.”  
That keeps the implementation scientifically falsifiable.  
**So my recommendation is: Jacobian/eigenvalue extraction → synthetic eigenvalue tests → stateful gate → SBLToolkit/GABLS3 integration → 24-h A/B validation.**  
If SBLToolkit.jl is already in your Library, I can inspect its existing GABLS3 driver and Jacobian/state interfaces and design SBLGating.jl against the actual API rather than inventing an integration layer.  

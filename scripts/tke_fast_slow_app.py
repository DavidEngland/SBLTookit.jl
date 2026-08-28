import streamlit as st
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import brentq
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Set Page Config
st.set_page_config(
    page_title="2D Fast-Slow TKE Dynamics",
    page_icon="🌪️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .reportview-container {
        background: #f0f2f6
    }
    .main-title {
        font-size: 38px;
        font-weight: 800;
        color: #1E3A8A;
        margin-bottom: 5px;
    }
    .sub-title {
        font-size: 18px;
        color: #4B5563;
        margin-bottom: 25px;
    }
    .section-header {
        font-size: 22px;
        font-weight: 600;
        color: #1E3A8A;
        border-bottom: 2px solid #E5E7EB;
        padding-bottom: 8px;
        margin-top: 20px;
        margin-bottom: 15px;
    }
    .highlight-card {
        background-color: #F3F4F6;
        border-left: 5px solid #3B82F6;
        padding: 15px;
        border-radius: 4px;
        margin-bottom: 15px;
    }
</style>
""", unsafe_allow_html=True)

# App Title
st.markdown('<div class="main-title">2D Fast-Slow Dynamics of Turbulent Kinetic Energy and Shear</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Interactive exploration of the Stable Boundary Layer (SBL) Dynamical System under Geometric Singular Perturbation Theory (GSPT)</div>', unsafe_allow_html=True)

# Sidebar for controls
st.sidebar.title("🎛️ Model Controls")

# Parameter Sliders
with st.sidebar.expander("🛠️ Model Parameters", expanded=True):
    epsilon = st.slider("Scale Separation (ε)", 0.005, 0.50, 0.05, step=0.005, help="Ratio of fast TKE relaxation rate to slow momentum advection rate.")
    l = st.slider("Mixing Length (l)", 1.0, 30.0, 10.0, step=0.5, help="Length scale representing the size of turbulent eddies.")
    N2 = st.slider("Brunt-Väisälä Frequency² (N²)", 0.0, 0.02, 0.005, step=0.0005, format="%.4f", help="Square of Brunt-Väisälä frequency, representing atmospheric stability.")
    phi = st.slider("Buoyancy Scaling Parameter (φ)", 0.1, 2.0, 1.0, step=0.1, help="Dimensionless buoyancy scaling parameter.")
    c_b = st.slider("Buoyant Sink Efficiency (c_b)", 0.0, 1.5, 0.5, step=0.05, help="Dimensionless saturation coefficient for buoyancy sink.")
    alpha = st.slider("Regularization Parameter (α)", 0.01, 0.5, 0.1, step=0.01, help="Regularization parameter in buoyancy saturation term.")
    G = st.slider("Geostrophic Forcing (G)", 0.01, 1.5, 0.5, step=0.01, help="External geostrophic forcing driving mean wind shear.")
    c_1 = st.slider("Shear Damping Coeff. (c_1)", 0.1, 5.0, 1.80, step=0.1, help="Turbulent feedback on mean wind shear.")
    c_2 = st.slider("Linear Decay Coeff. (c_2)", 0.0, 0.5, 0.10, step=0.01, help="Frictional/linear decay of mean shear.")
    delta = st.number_input("Regularization Floor (δ)", value=1e-6, format="%.1e", help="Velocity scale smoothing constant to ensure C¹ continuity.")

with st.sidebar.expander("🚀 Initial Conditions & Sim Settings", expanded=True):
    E_0 = st.slider("Initial TKE (E₀)", 0.0, 3.0, 0.1, step=0.05)
    S_0 = st.slider("Initial Shear (S₀)", 0.0, 1.5, 0.5, step=0.05)
    t_max = st.slider("Simulation Time (Seconds)", 10, 300, 100, step=10)
    solver_method = st.selectbox("Solver Method", ["Radau", "BDF", "RK45"], index=0, help="Radau or BDF are recommended as the fast-slow system is highly stiff.")

# Mathematical formulation definitions
def fast_subsystem_vector_field(E, S):
    # f(E, S) = l * sqrt(E + delta) * (S^2 - phi * N2 - c_b * N2 * E / (E + alpha)) - E^1.5 / l
    E_clamped = np.maximum(E, 0.0)
    production = S**2
    buoyancy = phi * N2 + c_b * N2 * E_clamped / (E_clamped + alpha)
    f = l * np.sqrt(E_clamped + delta) * (production - buoyancy) - (E_clamped**1.5) / l
    return f

def slow_subsystem_vector_field(E, S):
    # g(E, S) = G - c_1 * E * S - c_2 * S
    E_clamped = np.maximum(E, 0.0)
    return G - c_1 * E_clamped * S - c_2 * S

def local_fast_eigenvalue(E, S):
    # lambda_f = f_E(E, S)
    E_clamped = np.maximum(E, 1e-12)
    term1 = l * (S**2 - phi * N2) / (2.0 * np.sqrt(E_clamped + delta))
    term2 = c_b * N2 * l * (E_clamped / (2.0 * np.sqrt(E_clamped + delta) * (E_clamped + alpha)) + alpha * np.sqrt(E_clamped + delta) / (E_clamped + alpha)**2)
    term3 = 1.5 * np.sqrt(E_clamped) / l
    return term1 - term2 - term3

def ode_system(t, y):
    E, S = y
    # dy/dt = [f(E, S) / epsilon, g(E, S)]
    dE_dt = fast_subsystem_vector_field(E, S) / epsilon
    dS_dt = slow_subsystem_vector_field(E, S)
    return [dE_dt, dS_dt]

# Run simulation
t_eval = np.linspace(0, t_max, int(t_max * 10))
y0 = [E_0, S_0]
sol = integrate.solve_ivp(ode_system, (0, t_max), y0, method=solver_method, t_eval=t_eval)

E_trajectory = sol.y[0]
S_trajectory = sol.y[1]
t_trajectory = sol.t

# Find the unique fixed point
def find_fixed_point_helper(g_val):
    def eq(E):
        if E < 0:
            return -g_val / c_2 if c_2 > 0 else -1e3
        S_M0_val = np.sqrt(phi * N2 + c_b * N2 * E / (E + alpha) + (E**1.5) / (l**2 * np.sqrt(E + delta)))
        S_slow_val = g_val / (c_1 * E + c_2) if (c_1 * E + c_2) > 0 else 1e5
        return S_M0_val - S_slow_val
    
    # Check physical transition boundary: if G < c_2 * sqrt(phi * N2), fixed point is at E* = 0, S* = G / c_2
    S_crit = np.sqrt(phi * N2)
    if g_val <= c_2 * S_crit:
        return 0.0, g_val / c_2 if c_2 > 0 else 0.0
    
    try:
        E_max = 50.0
        while eq(E_max) < 0:
            E_max *= 5.0
        E_star = brentq(eq, 0.0, E_max)
        S_star = g_val / (c_1 * E_star + c_2)
        return E_star, S_star
    except Exception:
        # Fallback to analytical laminar state
        return 0.0, g_val / c_2 if c_2 > 0 else 0.0

E_star, S_star = find_fixed_point_helper(G)
G_crit = c_2 * np.sqrt(phi * N2)

# Create layout tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📈 Phase Space Portrait", 
    "🕒 Time Series Evolution", 
    "🛡️ Stability & Eigenvalues", 
    "🔀 Forcing Bifurcation",
    "📚 Theoretical Background"
])

# --- TAB 1: PHASE SPACE ---
with tab1:
    st.markdown('<div class="section-header">Phase Space Analysis</div>', unsafe_allow_html=True)
    st.write(
        "The phase space displays the coupled evolution of **Turbulent Kinetic Energy (E)** and **Wind Shear (S)**. "
        "It reveals how the system's trajectories rapidly relax onto the attracting critical manifold, "
        "then slide along it toward the steady state (fixed point)."
    )
    
    # Compute critical manifold (M_0)
    E_grid = np.linspace(0.0, np.maximum(np.max(E_trajectory) * 1.2, 2.5), 500)
    S_M0 = np.sqrt(phi * N2 + c_b * N2 * E_grid / (E_grid + alpha) + (E_grid**1.5) / (l**2 * np.sqrt(E_grid + delta)))
    
    # Compute slow nullcline
    S_slow = G / (c_1 * E_grid + c_2)
    
    # Plotly phase space figure
    fig_phase = go.Figure()
    
    # Add Critical Manifold M_0
    fig_phase.add_trace(go.Scatter(
        x=E_grid, y=S_M0,
        mode='lines',
        name='Critical Manifold M₀ (f=0)',
        line=dict(color='royalblue', width=3),
        hovertemplate='E: %{x:.4f}<br>S (M₀): %{y:.4f}<extra></extra>'
    ))
    
    # Add Slow Nullcline
    fig_phase.add_trace(go.Scatter(
        x=E_grid, y=S_slow,
        mode='lines',
        name='Slow Nullcline (g=0)',
        line=dict(color='orange', width=2, dash='dash'),
        hovertemplate='E: %{x:.4f}<br>S (Slow Null): %{y:.4f}<extra></extra>'
    ))
    
    # Add Trajectory
    fig_phase.add_trace(go.Scatter(
        x=E_trajectory, y=S_trajectory,
        mode='lines+markers',
        name='Integrated Trajectory',
        line=dict(color='firebrick', width=2),
        marker=dict(size=4, color='firebrick', symbol='circle'),
        hovertemplate='t: %{text:.2f}<br>E: %{x:.4f}<br>S: %{y:.4f}<extra></extra>',
        text=t_trajectory
    ))
    
    # Add Start Point
    fig_phase.add_trace(go.Scatter(
        x=[E_0], y=[S_0],
        mode='markers',
        name='Initial Condition (y₀)',
        marker=dict(color='green', size=12, symbol='star'),
        hovertemplate='Initial State:<br>E₀: %{x:.4f}<br>S₀: %{y:.4f}<extra></extra>'
    ))
    
    # Add Fixed Point
    fig_phase.add_trace(go.Scatter(
        x=[E_star], y=[S_star],
        mode='markers',
        name='Steady-State Fixed Point',
        marker=dict(color='black', size=12, symbol='x'),
        hovertemplate='Fixed Point:<br>E*: %{x:.4f}<br>S*: %{y:.4f}<extra></extra>'
    ))
    
    # Flow direction quiver/vector field (optional/toggle)
    show_flow = st.checkbox("Show Fast Vector Field (f(E, S) Flow)", value=True)
    if show_flow:
        E_vec_grid = np.linspace(0.0, np.maximum(np.max(E_trajectory) * 1.2, 2.5), 15)
        S_vec_grid = np.linspace(0.0, np.maximum(np.max(S_trajectory) * 1.2, 1.2), 15)
        EE, SS = np.meshgrid(E_vec_grid, S_vec_grid)
        
        # We only plot fast direction vectors to show the rapid relaxation
        f_val = fast_subsystem_vector_field(EE, SS) / epsilon
        
        # Normalize for visualization
        norm = np.abs(f_val)
        norm[norm == 0] = 1.0
        dx = f_val / norm * 0.05
        
        for i in range(len(E_vec_grid)):
            for j in range(len(S_vec_grid)):
                if np.isnan(dx[j, i]) or np.isinf(dx[j, i]):
                    continue
                fig_phase.add_trace(go.Scatter(
                    x=[EE[j, i], EE[j, i] + dx[j, i]],
                    y=[SS[j, i], SS[j, i]],
                    mode='lines',
                    line=dict(color='rgba(128, 128, 128, 0.3)', width=1),
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    fig_phase.update_layout(
        title="2D Phase Space Portrait (E vs S)",
        xaxis_title="Turbulent Kinetic Energy E (m² s⁻²)",
        yaxis_title="Mean Shear S (s⁻¹)",
        legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.8)"),
        margin=dict(l=40, r=40, t=50, b=40),
        height=600,
        hovermode="closest",
        template="plotly_white"
    )
    
    st.plotly_chart(fig_phase, use_container_width=True)
    
    # State values card
    st.markdown('<div class="highlight-card">', unsafe_allow_html=True)
    col1, col2, col3 = st.columns(3)
    col1.metric(label="Steady State TKE (E*)", value=f"{E_star:.4f} m² s⁻²")
    col2.metric(label="Steady State Shear (S*)", value=f"{S_star:.4f} s⁻¹")
    col3.metric(label="Critical Forcing Threshold (G_crit)", value=f"{G_crit:.4f} s⁻²")
    st.markdown('</div>', unsafe_allow_html=True)


# --- TAB 2: TIME SERIES ---
with tab2:
    st.markdown('<div class="section-header">Time Series Evolution</div>', unsafe_allow_html=True)
    st.write(
        "This view shows how **Turbulent Kinetic Energy (E)** and **Wind Shear (S)** evolve over time. "
        "Notice how the fast variable (TKE) exhibits a very steep, rapid transition initially, while "
        "the slow variable (Shear) changes gradually over a much longer timescale."
    )
    
    fig_time = make_subplots(specs=[[{"secondary_y": True}]])
    
    # TKE Trace
    fig_time.add_trace(
        go.Scatter(x=t_trajectory, y=E_trajectory, name="TKE E(t)", line=dict(color="firebrick", width=2.5)),
        secondary_y=False,
    )
    
    # Shear Trace
    fig_time.add_trace(
        go.Scatter(x=t_trajectory, y=S_trajectory, name="Shear S(t)", line=dict(color="navy", width=2.5, dash="dash")),
        secondary_y=True,
    )
    
    fig_time.update_layout(
        title_text="Coupled Fast-Slow Time Series Dynamics",
        xaxis_title="Time (s)",
        template="plotly_white",
        height=500,
        legend=dict(x=0.02, y=0.98)
    )
    
    fig_time.update_yaxes(title_text="<b>TKE E</b> (m² s⁻²)", color="firebrick", secondary_y=False)
    fig_time.update_yaxes(title_text="<b>Shear S</b> (s⁻¹)", color="navy", secondary_y=True)
    
    st.plotly_chart(fig_time, use_container_width=True)


# --- TAB 3: STABILITY ---
with tab3:
    st.markdown('<div class="section-header">Stability and Hyperbolicity Analysis</div>', unsafe_allow_html=True)
    st.write(
        "Geometric Singular Perturbation Theory dictates that the **eigenvalue of the fast subsystem (λ_f)** "
        "governs whether the system is attracted to or repelled by the critical manifold. "
        "If **λ_f < 0**, the manifold is **normally hyperbolic and attracting**."
    )
    
    # Calculate lambda_f along the integrated trajectory
    lam_trajectory = local_fast_eigenvalue(E_trajectory, S_trajectory)
    
    fig_stab = make_subplots(rows=1, cols=2, subplot_titles=(
        "Fast Eigenvalue λ_f(t) over Time",
        "Phase Trajectory Colored by Stability (λ_f)"
    ))
    
    # 1. Lambda_f time series
    fig_stab.add_trace(go.Scatter(
        x=t_trajectory, y=lam_trajectory,
        mode='lines',
        name='λ_f(t)',
        line=dict(color='darkviolet', width=2.5)
    ), row=1, col=1)
    # Add a horizontal line at 0 for threshold
    fig_stab.add_trace(go.Scatter(
        x=[0, t_max], y=[0, 0],
        mode='lines',
        name='Stability Threshold (λ_f = 0)',
        line=dict(color='black', width=1, dash='dot')
    ), row=1, col=1)
    
    # 2. Phase space colored by lambda_f
    fig_stab.add_trace(go.Scatter(
        x=E_trajectory, y=S_trajectory,
        mode='markers+lines',
        name='Trajectory (colored)',
        marker=dict(
            size=6,
            color=lam_trajectory,
            colorscale='Viridis',
            colorbar=dict(title="λ_f", x=1.02),
            showscale=True
        ),
        hovertemplate='E: %{x:.4f}<br>S: %{y:.4f}<br>λ_f: %{marker.color:.4e}<extra></extra>'
    ), row=1, col=2)
    
    # Add Critical Manifold for reference in phase plot
    fig_stab.add_trace(go.Scatter(
        x=E_grid, y=S_M0,
        mode='lines',
        name='M₀ (f=0)',
        line=dict(color='grey', width=1.5, dash='dash')
    ), row=1, col=2)
    
    fig_stab.update_layout(
        template="plotly_white",
        height=500,
        showlegend=False
    )
    
    fig_stab.update_xaxes(title_text="Time (s)", row=1, col=1)
    fig_stab.update_yaxes(title_text="Local Fast Eigenvalue λ_f", row=1, col=1)
    fig_stab.update_xaxes(title_text="TKE E (m² s⁻²)", row=1, col=2)
    fig_stab.update_yaxes(title_text="Mean Shear S (s⁻¹)", row=1, col=2)
    
    st.plotly_chart(fig_stab, use_container_width=True)
    
    # Stability text card
    st.markdown("""
    <div class="highlight-card">
    <h4>🔍 Interpretation of Stability</h4>
    <ul>
        <li><b>λ_f < 0 (Attracting):</b> The trajectory is exponentially drawn toward the manifold. Fast perturbations decay instantly, preserving the dimensional reconciliation.</li>
        <li><b>λ_f = 0 (Fold Boundary / Loss of Hyperbolicity):</b> This marks the boundary where the fast dynamics lose stability. Here, the system undergoes an explosive fast relaxation transition (collapse of active turbulence or rapid onset).</li>
        <li><b>λ_f > 0 (Repelling):</b> Trajectories are actively blown away from this region, preventing the system from maintaining quasi-equilibrium.</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)


# --- TAB 4: BIFURCATION ---
with tab4:
    st.markdown('<div class="section-header">Forcing Bifurcation Analysis</div>', unsafe_allow_html=True)
    st.write(
        "How does the steady-state of the atmosphere respond to external geostrophic forcing ($G$)? "
        "This diagram dynamically calculates the steady-state TKE ($E^*$) and Wind Shear ($S^*$) as a function of $G$, "
        "allowing you to see the critical threshold where the boundary layer transitions from laminar to turbulent."
    )
    
    # Generate G values array
    G_range = np.linspace(0.001, 1.2, 150)
    E_stars = []
    S_stars = []
    
    for g_val in G_range:
        e_st, s_st = find_fixed_point_helper(g_val)
        E_stars.append(e_st)
        S_stars.append(s_st)
        
    fig_bifur = make_subplots(rows=1, cols=2, subplot_titles=(
        "Steady-State TKE E* vs Geostrophic Forcing G",
        "Steady-State Shear S* vs Geostrophic Forcing G"
    ))
    
    # TKE vs G
    fig_bifur.add_trace(go.Scatter(
        x=G_range, y=E_stars,
        mode='lines',
        name='TKE E*',
        line=dict(color='firebrick', width=3)
    ), row=1, col=1)
    # Add vertical line at current G
    fig_bifur.add_trace(go.Scatter(
        x=[G, G], y=[0, max(E_stars)*1.1],
        mode='lines',
        name='Selected Forcing (G)',
        line=dict(color='grey', width=1.5, dash='dash')
    ), row=1, col=1)
    
    # Shear vs G
    fig_bifur.add_trace(go.Scatter(
        x=G_range, y=S_stars,
        mode='lines',
        name='Shear S*',
        line=dict(color='navy', width=3)
    ), row=1, col=2)
    # Add vertical line at current G
    fig_bifur.add_trace(go.Scatter(
        x=[G, G], y=[0, max(S_stars)*1.1],
        mode='lines',
        name='Selected Forcing (G)',
        line=dict(color='grey', width=1.5, dash='dash')
    ), row=1, col=2)
    
    fig_bifur.update_layout(
        template="plotly_white",
        height=500,
        showlegend=False
    )
    
    fig_bifur.update_xaxes(title_text="Geostrophic Forcing G (s⁻²)", row=1, col=1)
    fig_bifur.update_yaxes(title_text="Steady-State TKE E* (m² s⁻²)", row=1, col=1)
    fig_bifur.update_xaxes(title_text="Geostrophic Forcing G (s⁻²)", row=1, col=2)
    fig_bifur.update_yaxes(title_text="Steady-State Wind Shear S* (s⁻¹)", row=1, col=2)
    
    st.plotly_chart(fig_bifur, use_container_width=True)
    
    st.markdown(f"""
    <div class="highlight-card">
    <h4>📉 Critical Transition Dynamics</h4>
    <p>For your current parameter configuration, the critical geostrophic forcing threshold is <b>G_crit = {G_crit:.5f} s⁻²</b>.</p>
    <ul>
        <li><b>Weak Forcing (G ≤ G_crit):</b> The geostrophic forcing is too weak to sustain shear-driven turbulence. The stable boundary layer settles into a completely laminar state (<b>E* = 0</b>, <b>S* = G/c₂ = {G/c_2 if c_2 > 0 else 0:.4f} s⁻¹</b>).</li>
        <li><b>Strong Forcing (G > G_crit):</b> Shear-driven production overcomes the buoyant consumption. A transcritical-type bifurcation occurs, and a stable turbulent branch emerges (<b>E* > 0</b>), which acts as a powerful feedback mechanism keeping wind shear under control.</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)


# --- TAB 5: THEORETICAL BACKGROUND ---
with tab5:
    st.markdown('<div class="section-header">Theoretical Background and System Equations</div>', unsafe_allow_html=True)
    
    st.write(
        "This application models the coupled fast-slow dynamics of the Stable Boundary Layer (SBL) using "
        "Geometric Singular Perturbation Theory (GSPT). This specific formulation resolves historical "
        "dimensional and mathematical discrepancies in the atmospheric boundary layer by introducing the "
        "smoothing parameter $\\delta$ directly into the TKE budget vector field."
    )
    
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.markdown("### 🌀 Fast Subsystem (TKE)")
        st.latex(r"""
        \epsilon \frac{dE}{dt} = l\sqrt{E+\delta}\left(S^2 - \phi N^2 - c_b N^2 \frac{E}{E+\alpha}\right) - \frac{E^{3/2}}{l}
        """)
        st.markdown("""
        * **$E$**: Turbulent Kinetic Energy (TKE) - fast variable.
        * **$\\epsilon$**: Scale separation ratio ($0 < \\epsilon \\ll 1$).
        * **$l$**: Mixing length scale.
        * **$\\delta$**: Regularization floor ($10^{-6}$), smoothing the vector field at $E \\to 0$.
        * **$S^2$**: TKE production due to wind shear.
        * **$N^2$**: Buoyancy frequency (sink under stable stratification).
        * **$c_b \\frac{E}{E+\\alpha}$**: Buoyant sink saturation modeling TKE damping.
        """)
        
    with col_right:
        st.markdown("### 💨 Slow Subsystem (Shear)")
        st.latex(r"""
        \frac{dS}{dt} = G - c_1 E S - c_2 S
        """)
        st.markdown("""
        * **$S$**: Mean wind shear - slow variable.
        * **$G$**: External geostrophic wind forcing.
        * **$c_1 E S$**: Shear destruction due to turbulent momentum flux.
        * **$c_2 S$**: Frictional and linear momentum decay.
        """)
        
    st.markdown("### 📐 Exact Jacobian Derivative (Local Fast Eigenvalue)")
    st.write(
        "Differentiating the fast vector field $f(E, S)$ with respect to TKE ($E$) yields the exact local "
        "fast eigenvalue $\\lambda_f \\equiv f_E$, which governs the normal hyperbolicity of the system:"
    )
    st.latex(r"""
    \lambda_f = \frac{l\left(S^2 - \phi N^2\right)}{2\sqrt{E+\delta}} - c_b N^2 l \left[ \frac{E}{2\sqrt{E+\delta}(E+\alpha)} + \frac{\alpha\sqrt{E+\delta}}{(E+\alpha)^2} \right] - \frac{3\sqrt{E}}{2l}
    """)
    
    st.markdown("### 📌 Saddle-Node Fold Tracking via Newton-Raphson")
    st.write(
        "The coordinates of saddle-node fold points $(E_{\\text{fold}}, S_{\\text{fold}})$ are "
        "calculated by solving the 2D root-finding problem where both the fast field and its derivative vanish simultaneously:"
    )
    st.latex(r"""
    \mathbf{F}(E, S) = \begin{bmatrix} f(E, S) \\ f_E(E, S) \end{bmatrix} = \mathbf{0}
    """)
    st.write(
        "The system iterates using the multi-dimensional Newton-Raphson method with symbolic Jacobian updates "
        "until convergence is reached under $\\Vert\\mathbf{F}\\Vert_2 < 10^{-8}$:"
    )
    st.latex(r"""
    \mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} - \mathbf{J}_{\mathbf{F}}^{-1} \mathbf{F}(\mathbf{x}^{(k)})
    """)

# Instructions on running
st.sidebar.markdown("---")
st.sidebar.markdown("""
### 🏃 How to Run Locally:
1. Copy this code and save it as `app.py`
2. Install dependencies:
   ```bash
   pip install streamlit numpy scipy plotly
   ```
3. Run the app:
   ```bash
   streamlit run app.py
   ```
""")

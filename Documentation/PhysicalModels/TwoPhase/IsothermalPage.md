The isothermal two-fluid model
==============================

The model consists in the phasic mass and momentum balance equations.

The main unknowns are $\alpha$, $P$, $\vec{u}_g$, $\vec{u}_l$. The model uses stiffened gas laws $p_g(\rho_g)$ and  $p_l(\rho_l)$ valid around the saturation points $(P=1 bar, T=373K)$ or $(P=155 bars, T=618K)$ depending on the value of the enum \ref pressureEstimate.

The subscript $k$ stands for $l$ for the liquid phase and $g$ for the gas phase. The common
averaged pressure of the two phases is denoted by $p$. 

In our model, pressure equilibrium between the two phases is postulated, and the resulting system to solve is:
$$
\left\{
\begin{array}{ccll}
 \frac{\partial m_g}{\partial t}& +& \nabla \cdot \vec{q}_g &= 0,\\[1.5ex]
\frac{\partial m_l}{\partial t} &+ &\nabla \cdot \vec{q}_l &= 0,\\[1.5ex]
\frac{\partial \vec{q}_g}{\partial t}& +& \nabla \cdot (\vec{q}_g\otimes\frac{\vec{q}_g}{m_g})+ \alpha_g \vec\nabla p&\\[1.5ex] 
 &+&\Delta p \nabla \alpha_g -\nu_g\Delta \vec{u}_g &= m_g\vec{g}-K\rho_m||\vec{u}_g-\vec{u}_l||(\vec{u}_g-\vec{u}_l)-K_s\delta(x)m_g||\vec{u}_g||\vec{u}_g\\[1.5ex]
\frac{\partial \vec{q}_l}{\partial t}& +& \nabla \cdot (\vec{q}_l\otimes\frac{\vec{q}_l}{m_l})+ \alpha_l \vec\nabla p&\\[1.5ex]
&+&\Delta p \nabla \alpha_l -\nu_l\Delta \vec{u}_l &= m_l\vec{g}-K\rho_m||\vec{u}_l-\vec{u}_g||(\vec{u}_l-\vec{u}_g)-K_s\delta(x)m_l||\vec{u}_l||\vec{u}_l\\
\end{array}
\right. \label{eq: isentropictwophase}
$$

Here :
- $\nu_k$ is the viscosity of phase $k$, set by \ref IsothermalTwoFluid::setViscosity
- $\Delta p$ denotes the pressure default $p-p_k$ between the bulk average pressure and the interfacial average pressure.
- $\vec g$ the gravity vector (\ref IsothermalTwoFluid::setGravity)
- $K$ the interphase friction coefficient (\ref IsothermalTwoFluid::setDragCoeffs),
- $K_s(\vec x)$ the singular friction function, $\delta_s(\vec x)$ the Dirac delta function with support on the set $s$ (\ref IsothermalTwoFluid::setPressureLossField).

where $ 
 \left\{\begin{array}{clc}
	\alpha_g +\alpha_l &=& 1 \\[1.5ex]
	m_k &=& \alpha_k \rho_k \\[1.5ex]
	\vec{q}_k &=& \alpha_k \rho_k \vec{u}_k \\[1.5ex]
        
	\end{array}\right..
$


* The class : \ref IsothermalTwoFluid implements the isentropic two-fluid model

* \subpage ExampleIsothermalPage "Here are C and Python example scripts using the isothermal two-fluid model"	


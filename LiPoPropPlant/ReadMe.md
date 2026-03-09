# GOAL
`an effect of battery discharge simulation`

# Math model (idea)
```
Battery:
	Input:
		- required current (I_bat)
	Output:
		- Voltage (Vt)
		- Current (I_bat)
	Internal state:
		- State-of-Charge
	Description:
		by system of ODEs
			d(SoC)/dt = -I^k / Cp
			dVp/dt = -Vp / (Rp*Cpol) + I/Cp
			Vt = Voc(SoC) - I * R0 - Vp
			where
			Cp is Peukert Capacity
			Cpol is Polarization Capacitance
			k is he Peukert constant (typically 1.05 to 1.3 depending on battery chemistry).
			R0 is Internal ohmic resistance.
ESC:
	Input: 
		- Voltage (Vt)
		- Request motor current (I_motor)
		- throttle signal (delta) from flight computer
		- motor current (I_motor) from motor
	Output:
		- V_eff (effective voltage to motor)
		- battery current (I_bat) to motor and to battery as a feedback
	Internal state:
		- I_qiescent is the tiny current used to power the ESC's own microchip
		- Resc is internal resistance
	Description:
		by algebraic equations
			Veff = delta * Vt - I_motor * Resc
			I_bat = delta * I_motor + I_qiescent

Motor:
	Input:
		- effective voltage Veff
		- T_load
		- current (I_bat)
	Output:
		- motor RPM 
		- required current to produce the RPM (I_motor)
	Internal state:
		- Ke is Back-EMF Constant
		- Rm is represents the resistance of the copper windings
	Description:
		by algebraic equations
			I_motor = (Veff - Ke*omega) / Rm (from Veff = Ke*omega + I_bat * Rm)
			T_motor = Kt * I_motor
			where omega = pi * RPM / 30

Fan:
	Input:
		- omega
	Output:
		- thrust, T_load is a breaking force
	Internal state:
		- number of blades
		- diameter
		- aerodynamic characteristics
	Description:
		by algebraic equations
			T -> Ct * omega^2
			J*omega_dot = T_motor - T_load - T_friction
			where J is the moment of inertia of the rotating assembly

```
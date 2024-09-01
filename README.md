# General MPC builder with CasADi

## Overview

The routine _mpc_build_ is written to build a generic MPC formulation with CasADi. The input arguments are all (or, at least, most of) the possible configurations of an MPC problem. All the elements of such a configuration are to be passed under a structure of options. 

## Options

### Modeling

The following lists the options possible (some are mandatory):

* Number of states, inputs and the prediction horizon: _opt.n_states_, _opt.n_controls_ and _opt.N_, respectively.
* Nature of the model: _opt.model.type_, which should be either 'linear', or 'nonlinear'.
  * If _opt.type = linear_, the defining matrices A and B should be provided through _opt.model.A_ and _opt.model.B_.
  * If _opt.type = nonlinear_, the respective function, being a function handle, must be provided through _opt.model.fun_. Note that the independent variables must ordered as @(x,u).
    * If the model described in _opt.model.function_ is *continuous*, a field _opt.continuous_model_ is expected. This field contains_opt.continuous_model.integration_, which defines the integration scheme to be undertaken to cast the prediction and the options are "euler" and "RK4" (forth-order Runge Kutta), and _opt.dt_, which is the time step.

Example:

```matlab
opt.n_states   	= 2;
opt.n_controls 	= 2;
opt.model.type 	= 'nonlinear'
opt.model.fun  	= @(x,u) [x(1)^2+u(1);x(2)^2+u(2)];
opt.continuous_model.integration = 'euler';
opt.dt		= 0.1;
```

### Constraints

These terms relate to general constraints to be imposed to the optimization problem. All kinds of constraints are gathered in _opt.constraints_. The types of constraints supported are:

* Varible-wise bounds, provided through _opt.constraints.states.upper_ and __opt.constraints.states.lower__ for states, and _opt.constraints.control.upper_ and __opt.constraints.control.lower for the control inputs. These bounds must be given as a vector of appropriate dimensions (_i.e._, _opt.n_states_ and _opt.n_controls_).
* The state constraints can be polyhedral, _i.e._, *Ax\leq , and such an argument is to be provided through _opt.constraints.polyhedral_. This argument is expected to be composed of matrices $A$ and $b$.
  * Note that _Polyhedron_ objects, as those created by the MPT toolbox, are acceptable.

Terminal constraints (polyhedral or end-point) are possible, and should be provided through _opt.constraints.terminal_:

* A set, in the form of a polyhedron (as above), is expected in _opt.constraints.terminal.set_, and the respective matrice A and b are needed (_i.e._, $Ax(N)\leq b$).
* If any other parameters (seem as decision variables for the optimization problem) are to be terminally-constrained, they mst be listed in _opt.constraints.terminal.parameters_.

If any other parameter/decision variable is to be constrained (not terminally), there is an option _opt.constraints.parameters.variables_. A list of all constrained variables are expected. The corresponding bounds of such variables are expected in _opt.constraints.parameters.upper_ and _opt.constraints.parameters.lower_.

### Stage costs

These terms relate to state and control variables at each step over the prediction horizon (not at end-point, see *Terminal elements* below). The stage cost components are to be provided through _opt.costs.stage_:

* The stage cost function is provided through _opt.costs.stage.function_, which takes a function handle as arguments. This handle takes arguments @(x,u,extra), where "extra" is mandatory even if no other variables is taken into account.
  * If no _opt.costs.stage.function_ is provided, the simple linear-quadratic stage cost is considered. This one requires the weighting matrices _opt.costs.stage.Q_ and _opt.costs.stage.R_ as numerical matrices.
  * If there are no extra parameters to be considered in the stage cost function, simply add `0*extra`.
* If any other parameters (supposedly decison variables for the optimization problem) are considered in the stage cost function, it should be listed in the field _opt.costs.stage.parameters_.

### Terminal costs

These terms relate to state and control variables at the end of the prediction horizon (terminal point). The terminal ingredients are to be provided through _opt.costs.terminal_:

* The terminal cost function is provided through _opt.costs.terminal.function_, which takes a function handle as argument. This handle must take arguments @(x,u,extra), where "extra" are *optional*.
  *  _opt.costs.terminal.function_ is not mandatory.
* If any parameters (other than state $x$) are considered in the terminal cost function, it must be listed in the filed _opt.costs.stage.parameters_

### General constraints

These terms relate to general state and control constraints, which can be, for instance, a nonlinear function of the statem, or dependent on external parameters (therefore, time-varying).

### User inputs
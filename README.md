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

### Stage elements (costs and constraints)

These terms relate to state and control variables at each step over the prediction horizon (not at end-point, see *Terminal elements* below). The stage ingredients are to be provided through _opt.costs.stage_:

* The stage cost function is provided through _opt.costs.stage.function_, which takes a function handle as arguments. This handle takes arguments @(x,u,extra), where "extra" is mandatory even if no other variables is taken into account.
  * If no _opt.costs.stage.function_ is provided, the simple linear-quadratic stage cost is considered. This one requires the weighting matrices _opt.costs.stage.Q_ and _opt.costs.stage.R_ as numerical matrices.
  * If there are no extra parameters to be considered in the stage cost function, simply add `0*extra`.
* If any other parameters (supposedly decison variables for the optimization problem) are considered in the stage cost function, it should be listed in the field _opt.costs.stage.parameters_.

### Terminal elements (costs and constraints)

These terms relate to state and control variables at the end of the prediction horizon (terminal point). The terminal ingredients are to be provided through _opt.costs.terminal_:

### General constraints

### User inputs
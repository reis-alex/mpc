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

```
opt.n_states   	= 2;
opt.n_controls 	= 2;
opt.model.type 	= 'nonlinear'
opt.model.fun  	= @(x,u) [x(1)^2+u(1);x(2)^2+u(2)];
opt.continuous_model.integration = 'euler';
opt.dt		= 0.1;
```

### Stage elements (costs and constraints)

### Terminal elements (costs and constraints)

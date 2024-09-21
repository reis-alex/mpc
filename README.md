# General MPC builder with CasADi

## Overview

The routine _mpc_build_ is written to build a generic MPC formulation with CasADi. The input arguments are all (or, at least, most of) the possible configurations of an MPC problem. All the elements of such a configuration are to be passed under a structure of options. 

1. [Modeling](#Modeling)
2. Constraints
	1. [General Constraints](#Constraints)
	2. [General Constraints](#General-constraints)
3. Costs
	1. [Stage Costs](#Stage-costs)
	2. [Terminal Costs](#Terminal-costs)
4. [User Inputs](#User-inputs)

Some examples are available.

1. [MPC for tracking linear systems using artificial references](https://github.com/reis-alex/mpc/blob/main/Examples/Example_Tracking.md):

## Options

### Modeling

The following lists the options possible (some are mandatory):

* Number of states, inputs and the prediction horizon: _opt.n_states_, _opt.n_controls_ and _opt.N_, respectively.
* Nature of the model: _opt.model.type_, which should be either 'linear', or 'nonlinear'.
  * If _opt.type = linear_, the defining matrices A and B should be provided through _opt.model.A_ and _opt.model.B_.
  * If _opt.type = nonlinear_, the respective function, being a function handle, must be provided through _opt.model.fun_. Note that the independent variables must ordered as @(x,u).
    * If the model described in _opt.model.function_ is *continuous*, a field _opt.continuous_model_ is expected. This field contains_opt.continuous_model.integration_, which defines the integration scheme to be undertaken to cast the prediction and the options are "euler" and "RK4" (fourth-order Runge Kutta), and _opt.dt_, which is the time step.

Example:

Define the system
$\dot{x}_1 = x_1^2+u_1, \quad \dot{x}_2 = x_2^2+u_2$
and discretize it, using first-order Euler, with a time step $dt=0.1$.

```matlab
opt.n_states   	= 2;
opt.n_controls 	= 2;
opt.model.type 	= 'nonlinear'
opt.model.fun  	= @(x,u) [x(1)^2+u(1);x(2)^2+u(2)];
opt.continuous_model.integration = 'euler';
opt.dt		= 0.1;
```

### Parameters

Parameters are any decision variable to the optimization problem, or any input that is to be provided to the MPC online (references, for instance). The structure ``opt.parameters`` gathers all parameters to be declared by name and by dimension, _e.g._,

```matlab
opt.parameters.name = {'xs','us','ref','vector','matrix'};
opt.parameters.dim = [2 1; 1 1; 2 1; 10 1; 10 2];
```

In the example above, we declare four variables: $x_s\in\mathbb{R}^2$,  $u_s\in\mathbb{R}$, _ref_$\in\mathbb{R}^2$, _vector_ $x_s\in\mathbb{R}^10$, and _matrix_ $x_s\in\mathbb{R}^{10\times 2}$. Note that each line in ``opt.parameters.dim`` gathers the dimensions for each variable *in the order* declared in ``opt.parameters.name``.

### Constraints

These terms relate to general constraints to be imposed to the optimization problem. All kinds of constraints are gathered in _opt.constraints_. The types of constraints supported are:

* Varible-wise bounds, provided through _opt.constraints.states.upper_ and _opt.constraints.states.lower_ for states, and _opt.constraints.control.upper_ and _opt.constraints.control.lower_ for the control inputs. These bounds must be given as a vector of appropriate dimensions (_i.e._, _opt.n_states_ and _opt.n_controls_).

Example: Let $x\in\mathbb{R}^2$ and $u\in\mathbb{R}$. Define box-like constraints $\vert x_k\vert \leq 1$ and $\vert u_k\vert \leq 0.1$:

```matlab
opt.constraints.states.upper =  [1 1]; 
opt.constraints.states.lower = -[1 1]; 
opt.constraints.control.upper =  0.1;
opt.constraints.control.lower = -0.1
```

* The state constraints can be polyhedral, _i.e._, $Ax\leq b$, and such an argument is to be provided through _opt.constraints.polyhedral_. This argument is expected to be composed of matrices $A$ and $b$.
  * Note that _Polyhedron_ objects, as those created by the [MPT3 toolbox](https://www.mpt3.org/), are acceptable.

Example: same as above, but with a polyhedral definition:

```matlab
A = vertcat(eye(opt.n_states),-eye(opt.n_states));
b = ones(opt.n_states*2,1);
opt.constraints.terminal.set.A = A;
opt.constraints.terminal.set.b = b;
```
or, equivalently, if one uses the [MPT3 toolbox](https://www.mpt3.org/) (this might simplifying plotting the constraint set later):
```matlab
X = Polyhedron('A',vertcat(eye(opt.n_states),-eye(opt.n_states)),'b',ones(opt.n_states*2,1));
opt.constraints.terminal.set = X;
```

Terminal constraints (polyhedral or end-point) are possible, and should be provided through _opt.constraints.terminal_:

* A set, in the form of a polyhedron (as above), is expected in _opt.constraints.terminal.set_, and the respective matrice A and b are needed (_i.e._, $Ax(N)\leq b$).
* If any other parameters (seem as decision variables for the optimization problem) are to be terminally-constrained, they mst be listed in _opt.constraints.terminal.parameters_.

If any other parameter/decision variable is to be constrained (not terminally), there is an option _opt.constraints.parameters.variables_. A list of all constrained variables are expected. The corresponding bounds of such variables are expected in _opt.constraints.parameters.upper_ and _opt.constraints.parameters.lower_.

### General constraints

These options relate to general state and control constraints, which can be, for instance, a nonlinear function of the state, or dependent on external parameters (therefore, time-varying).

*Example* (non-linear constraint): suppose $x\in\mathbb{R}^2$ and define a nonlinear constraint such as $x_1^2 + x_2^2 \leq 1$:

```matlab
opt.constraints.general = @(x) x(1)^2+x(2)^2-1;
```

*Example* (constraint dependent on external parameters): suppose $x\in\mathbb{R}^2$ (```opt.n_states=2```) and that the state is to be constrained by a polyhedral set that is updated with time, _e.g._, $A(k)x\leq b(k)$. 

First, since the inputs ($A$ and $b$) are exclusively used in the general constraint, one should declare it as such. Clearly, since the MPC solver will not be updated afterwards, these input parameters have *fixed dimensions*, which we consider $10$ in this example:

```matlab
n_rows = 10;
opt.input.general_constraints.vector.dim = n_rows;
opt.input.general_constraints.matrix.dim = [n_rows opt.n_states];
```

Now, declare the constraining function:

```matlab
opt.constraints.general.function = @(x,varargin) varargin{1}*x([1 3]) - varargin{2};
```

Finally, when passing it as arguments to the MPC solver, one does as follows:

```matlab
A = ones(10,2); % any matrix with proper dimensions
b = ones(10,1); % any vector with proper dimensions
args.p                  = [initial conditions; other_parameters; reshape(A,n_rows*opt.n_states,1); b];
```

One should not that the _reshape_ function above is necessary since CasADi only takes vectors as inputs, therefore it is needed to concatenate the columns of A into a single vector. If the dimensions of the $A$ and $b$ changes at each iteration, an option is to declare _n_rows_ as a higher value and, online, complete the remaining rows of $A$ and $b$ with zeros.

### Stage costs

These terms relate to state and control variables at each step over the prediction horizon (not at end-point, see *Terminal elements* below). The stage cost components are to be provided through _opt.costs.stage_:

* The stage cost function is provided through _opt.costs.stage.function_, which takes a function handle as arguments. This handle takes arguments @(x,u,extra), where "extra" is mandatory even if no other variables is taken into account.
  * If no _opt.costs.stage.function_ is provided, the simple linear-quadratic stage cost, _i.e._, $V(x_k) = \sum_{k=0}^{N-1} x_k^\top Q x_k + u_k^\top R u_k$, is considered. This one requires the weighting matrices _opt.costs.stage.Q_ and _opt.costs.stage.R_ as numerical matrices of proper dimensions.
  * If there are no extra parameters to be considered in the stage cost function, simply add `0*extra`.
* If any other parameters (supposedly decision variables for the optimization problem) are considered in the stage cost function, it should be listed in the field _opt.costs.stage.parameters_.

Example: consider the classical linear-quadratic cost for tracking a constant reference $p$ with a repulsion term regarding $\sigma$:

$$
\begin{equation*}
V(x_k,u_k) = \sum_{i=1}^{N-1} (x_k-p_k)^\top Q (x_k-p_k) + u_k^\top R u_k + 100*(x_k-\sigma)^2
\end{equation*}
$$
```matlab
p = [10;10];
sigma = [5;5]; 
Q = eye(opt.n_state);
R = eye(opt.n_controls);
opt.costs.stage.function = @(x,u,extra) (x-p)'*Q*(x-p) + u'*R*u + 100*(x_k-\sigma)^2 + extra*0;
```

### Terminal costs

These terms relate to state and control variables at the end of the prediction horizon (terminal point). The terminal ingredients are to be provided through _opt.costs.terminal_:

* The terminal cost function is provided through _opt.costs.terminal.function_, which takes a function handle as argument. This handle must take arguments @(x,u,extra), where "extra" are *optional*.
  *  _opt.costs.terminal.function_ is not mandatory.
* If any parameters (other than state $x$) are considered in the terminal cost function, it must be listed in the filed _opt.costs.stage.parameters_

Example: consider a terminal constraint to be $x_N\in\Omega$, where $\Omega$ is a polyhedral set described such as $\Omega.A*x(N)\leq \Omega.b$.

```matlab
opt.constraints.terminal.set.A = Omega.A;
opt.constraints.terminal.set.b = Omega.b;
```
Note that it can be similarly done using a Polyhedron object (see Section _Constraints_ above).


### User inputs

One can pass external inputs to the optimization problem through the option _opt.input.vector_. The only input which is coded by default is the initializing state for the prediction (the feedback value measured from the state). Note that the variables provided through this option *are not* decision variables, but should still be parsed as CasADi symbolic variables. 

Example: consider a linear quadratic cost for tracking a given constant reference:

```matlab
opt.costs.stage.function = @(x,u,extra) (x-extra)'*Q*(x-extra) + u'*R*u;
ref = SX.sym('ref',opt.n_states);
opt.costs.stage.parameters = ref;
opt.input.vector = ref;
```
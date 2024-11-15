# General MPC builder with CasADi

## Overview

The routine _mpc_build_ is written to build a generic MPC formulation with CasADi. The input arguments are all (or, at least, most of) the possible configurations of an MPC problem. All the elements of such a configuration are to be passed under a structure of options. 

1. [Modeling](#Modeling)
2. Constraints
	1. [General Constraints](#Constraints)
	2. [General Constraints](#General-constraints)
3. [Parameters](#Parameters)
4. Costs
	1. [Stage Costs](#Stage-costs)
	2. [Terminal Costs](#Terminal-costs)


Some examples are available.

1. [MPC for tracking linear systems using artificial references](https://github.com/reis-alex/mpc/blob/main/Examples/Example_Tracking.md):
2. [MPC and robotics: using an URDF to generate a prediction model](https://github.com/reis-alex/mpc/blob/main/Examples/Example_URDF.md):

## Options

All options for building the MPC controller need to be passed through the structure ```opt```, whose fields are described in the following.

### Modeling

* ```opt.n_states```: number of states
* ```opt.n_controls```: number of controlled inputs
* ```opt.model.function```: system model. There are two ways to define this function:
    	1. As a function handle, _e.g._, ``@(x,u) A*x+B*u`` (note that the order of $x$ and $u$ must be respected)
	2. As a CasADi function, depending on variables SX.sym. Note that in this case, the state and input vector must be provided through _opt.model.states_ and _opt.model.controls_.
* ```opt.continuous_model_ ```: if the mode is in continuous time, the field ```opt.continuous_model.integration``` must be provided with the integration scheme to be undertaken when casting the prediction. The options are "euler" and "RK4" (fourth-order Runge Kutta)
* ```opt.dt```: time step for the integration scheme

Example:

Define the system
$\dot{x}_1 = x_1^2+u_1, \quad \dot{x}_2 = x_2^2+u_2$
and discretize it, using first-order Euler, with a time step $dt=0.1$.

```matlab
opt.n_states   	= 2;
opt.n_controls 	= 2;
opt.model.function  	= @(x,u) [x(1)^2+u(1);x(2)^2+u(2)];
opt.continuous_model.integration = 'euler';
opt.dt		= 0.1;
```

The same example but using CasADi variables:

```matlab
opt.n_states   	= 2;
opt.n_controls 	= 2;
x1  = SX.sym('x1');
x2 = SX.sym('x2');
u1  = SX.sym('u1');
u2 = SX.sym('u2');
states = [x1;x2];
controls = [u1,u2];
opt.model.function  	= [x1^2+u1;x2^2+u2];
opt.continuous_model.integration = 'euler';
opt.dt		= 0.1;
```

### Parameters

Parameters are any decision variable to the optimization problem, or any input that is to be provided to the MPC online (references, for instance). The structure ``opt.parameters`` gathers all parameters to be declared by name and by dimension, _e.g._,

* ```opt.parameters.name```: holds the name of all parameters, independently of where they will appear in the MPC structure.
* ```opt.parameters.dim```: holds the dimension of each parameter, in the order they apppear in ```opt.parameters.name```. If the parameter is a vector, the dimension expected is ```[rows 1]```, whereas if it is a matrix, ```[rows columns]```.

Example: let the parameters be four decision variables: $x_s\in\mathbb{R}^2$,  $u_s\in\mathbb{R}$, _ref_$\in\mathbb{R}^2$, $b\in\mathbb{R}^10$, and $A\in\mathbb{R}^{10\times 2}$.

```matlab
opt.parameters.name = {'xs','us','ref','b','A'};
opt.parameters.dim = [2 1; 1 1; 2 1; 10 1; 10 2];
```

Later, to introduce these variables in the MPC (_e.g., in constraints or cost functions), one will just refer to the names listed in ``opt.parameters.name``. *Important:* if any parameter is to be used in constraints, they should be listed in the field ```opt.constraints.parameters.variables``` (see below).

### Constraints

The field _opt.constraints_ gathers the constraints to be imposed to the optimization problem.

* ```opt.constraints.state.upper``` and ```opt.constraints.state.lower```: variable-wise, upper and lower bounds for the states. The expected argument are vectors with the dimensions ```opt.n_states```.
* ```opt.constraints.controls.upper``` and ```opt.constraints.controls.lower```: variable-wise, upper and lower bounds for the control inputs. The expected argument are vectors with the dimensions ```opt.n_controls```.

Example: Let $x\in\mathbb{R}^2$ and $u\in\mathbb{R}$. Define box-like constraints $\vert x_k\vert \leq 1$ and $\vert u_k\vert \leq 0.1$:

```matlab
opt.constraints.states.upper =  [1 1]; 
opt.constraints.states.lower = -[1 1]; 
opt.constraints.control.upper =  0.1;
opt.constraints.control.lower = -0.1
```

* ```opt.constraints.polyhedral```: define polyhedral constraints _i.e._, $Ax\leq b$. This argument is expected to be composed of matrices $A$ and $b$.
  * Note that _Polyhedron_ objects, as those created by the [MPT3 toolbox](https://www.mpt3.org/), are acceptable.

Example: same as above, but with a polyhedral definition:

```matlab
A = vertcat(eye(opt.n_states),-eye(opt.n_states));
b = ones(opt.n_states*2,1);
opt.constraints.polyhedral.set.A = A;
opt.constraints.polyhedral.set.b = b;
```
or, equivalently, if one uses the [MPT3 toolbox](https://www.mpt3.org/) (this might simplifying plotting the constraint set later):
```matlab
X = Polyhedron('A',vertcat(eye(opt.n_states),-eye(opt.n_states)),'b',ones(opt.n_states*2,1));
opt.constraints.polyhedral.set = X;
```

* ```opt.constraints.parameters.variables```: gathers all parameters/decision variables that appear, somehow, in constraints.

* ```opt.constraints.terminal```: define terminal constraints (either polyhedral or end-point). The expected arguments are:
	* polyhedron (as above), passed directly through the field ```opt.constraints.terminal.set```
	* matrices $A$ and $b$, to be passed through the fields ```opt.constraints.terminal.set.A``` and ```opt.constraints.terminal.set.b```

* opt.constraints.terminal.parameters```: used if any other parameter/decision variable is to be used in the terminal constraint.

Example: consider a terminal constraint to be $x_N\in\Omega$, where $\Omega$ is a polyhedral set described such as $\Omega.A*x(N)\leq \Omega.b$.

```matlab
opt.constraints.terminal.set.A = Omega.A;
opt.constraints.terminal.set.b = Omega.b;
```
Note that it can be similarly done using a Polyhedron object (see Section _Constraints_ above).


### General constraints

One can also define (multiple) general constraints, for instance, as nonlinear functions of the state, or dependent on external parameters (therefore, time-varying). The corresponding fields are:

* ```opt.constraints.general.parameters```: list of all parameters used in the general constraints, if any.
* ```opt.constraints.general.function{i}```: holds the function describing the constraint. 
* ```opt.constraints.general.elements{i}```: indicate which elements of the predicted states are taken into account in the general constraint. It can be either ``` 'end' ```, for only the last predicted state, or ``` 'N' ```, for all the instances of the predicted states.
* ```opt.constraints.general.type{i}```: indicate if the constraint is an ``` 'equality' ``` or an ``` 'inequality' ```.

*Important remark:* note that some of the fields above *are arrays*, therefore ```{i}``` should be a integer number describing cardinality, starting from 1.

*Example* (non-linear constraint): suppose $x\in\mathbb{R}^2$ and define a nonlinear constraint such as $x_1^2 + x_2^2 - p \leq 1$:

```matlab
opt.constraints.general.parameters  = {'p'}
opt.constraints.general.function{1} = @(x,varargin) x(1)^2+x(2)^2-varargin{1}-1;
opt.constraints.general.elements{1} = 'N';
opt.constraints.general.type{1} = 'equality';
```

### Stage costs

These terms relate to state and control variables at each step over the prediction horizon (not at end-point, see *Terminal elements* below). The stage cost components are to be provided through _opt.costs.stage_:

* The stage cost function is provided through _opt.costs.stage.function_, which takes a function handle as arguments. This handle takes arguments @(x,u,varargin).
* If any other parameters (supposedly decision variables for the optimization problem) are considered in the stage cost function, it should be listed in the field _opt.costs.stage.parameters_.

Example: consider the classical linear-quadratic cost for tracking a user-input reference $p$ with a repulsion term regarding $\sigma$:

$$
\begin{equation*}
V(x_k,u_k) = \sum_{i=1}^{N-1} (x_k-p_k)^\top Q (x_k-p_k) + u_k^\top R u_k + 100*(x_k-\sigma)^2
\end{equation*}
$$
```matlab
p = [10;10];
sigma = [5;5]; 
opt.parameters.name = {'p'};
opt.parameters.dim = [opt.n_states 1];
Q = eye(opt.n_state);
R = eye(opt.n_controls);
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:})'*Q*(x-varargin{:}) + u'*R*u + 100*(x-\sigma)^2 ;
opt.costs.stage.parameters = {'p'};
```

### Terminal costs

These terms relate to state and control variables at the end of the prediction horizon (terminal point). The terminal ingredients are to be provided through _opt.costs.terminal_:

* The terminal cost function is provided through _opt.costs.terminal.function_, which takes a function handle as argument. This handle must take arguments @(x,varargin), where "varargin" gathers optional parameters.
* If _opt.costs.terminal.function_ is not given, the parser will simply constrain $x(N)\in\mathbb{X}_c$, being $\mathcal{X}_c$ the state constraint set for the whole prediction.
* If any parameters (other than state $x$) are considered in the terminal cost function, it must be listed in the filed _opt.parameters_ and then declared to _opt.costs.terminal.parameters_.

Example: consider the problem of steering the system to the closest steady-state $x_s$ approaching a given reference $ref$. To do so, one might weight the terminal point of the prediction plus a deviation error:
```matlab
P = ...;
T = 1000*P;
opt.parameters.name = {'xs','Ref'};
opt.parameters.dim = [opt.n_states 1; opt.n_states 1];
opt.costs.terminal.parameters = {'xs','Ref'};
opt.costs.terminal.function = @(x,varargin) (x-varargin{:}(1:4))'*P*(x-varargin{:}(1:4)) + (varargin{:}(1:4)-varargin{:}(5:8))'*T*(varargin{:}(1:4)-varargin{:}(5:8))
```

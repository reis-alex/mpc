function [usol,xpred,rest,cost,solx] = solve_mpc(opt,solver,args)

if length(args.p)~=length(args.vars{2})
    error('Wrong number of inputs (p).')
end
sol             = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
solx            = full(sol.x);
usol = reshape(solx(opt.n_states*(opt.N+1)+1:opt.n_states*(opt.N+1)+opt.N*opt.n_controls),opt.n_controls,opt.N);
xpred = reshape(solx(1:opt.N*opt.n_states),opt.n_states,opt.N);
rest = solx(opt.n_states*(opt.N+1)+opt.N*opt.n_controls+1:end);
cost = full(sol.f);
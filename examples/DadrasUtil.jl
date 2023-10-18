@inline function eom_dadras!(du,u,p,t)
    a = p[1];b=p[2];c=p[3]
    du[1] = a*u[1] - u[2]*u[3]+u[4]
    du[2] = u[1]*u[3]-b*u[2]
    du[3] = u[1]*u[2]-c*u[3]+u[1]*u[4]
    du[4] = -u[2]

    return nothing
end

function eom_dadras(u,p,t)
    a = p[1];b=p[2];c=p[3]
    du = zeros(4)
    du[1] = a*u[1] - u[2]*u[3]+u[4]
    du[2] = u[1]*u[3]-b*u[2]
    du[3] = u[1]*u[2]-c*u[3]+u[1]*u[4]
    du[4] = -u[2]
    return du
end

@inline function jac_dadras(J,u,p,t)
    a = p[1];b=p[2];c=p[3]
    J .= [a -u[3] -u[2] 1;
         u[3] -b u[1] 0;
         u[2]+u[4] u[1] -c u[1];
         0 -1 0 0 ]
end

@inline function rescale(x)
    return x/sqrt(norm(x))
end

@inline function rescale_inv(y)
    return y*norm(y)
end

"""
    function rescaled_eom(du,u,p,t)

Derivative of f(rescale_inv)
"""
@inline function rescaled_eom(u,p,t)
    v = rescale_inv(u)
    dv = eom_dadras(v,p,t)

    return ForwardDiff.jacobian(rescale, v) * dv
end

@inline function rescaled_eom!(du,u,p,t)
    v = rescale_inv(u)
    eom_dadras!(du,v,p,t)
    du .= ForwardDiff.jacobian(rescale, v) * du

    return nothing
end

function dadrasTimeSeries()
    p = [8;40;14.9]
    ic = [10.;1.;10.;1.]
    dt = .01
    ds = ContinuousDynamicalSystem(eom_dadras!, ic, p)
    sol,_ = trajectory(ds, 100000, Δt=dt)
    dat = Matrix(sol)'[:,:]
    m = 1000000
    resample_boxsize = .8 # NOTE: changed from 0.1
    resample_sb_r = .2
    Y_res, t_vec = resampleToDistance(ds, dt, dat[:,1:m], resample_boxsize; pp=rescale, sb_r=resample_sb_r, sb_fct=x->rescaled_eom(rescale(x), p, 0), max_depth=512, verbose=false)
    Y_res_rescaled = mapslices(rescale, Y_res, dims=[1])
    Z_res_rescaled = mapslices(x -> normalize(rescaled_eom(x, p, 0),2), Y_res_rescaled, dims=[1])

    return Y_res_rescaled, Z_res_rescaled, t_vec
end
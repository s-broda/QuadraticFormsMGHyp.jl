# computes the pdf of the 2SLS estimator under multivariate t errors
Random.seed!(0)
order = 2
n = 25
beta = 0.
k = 1
Z = randn(n, k)
In = diagm(ones(n))
Z = In[:, 1:k] * (Z'*Z)^(.5) # distribution depends on Z only through Z'Z if R=I
alpha = 0.
if abs(alpha) >= 1.
     b = 0.
else
     b = 1 / sqrt(1-alpha^2)
end
R = Toeplitz(alpha .^(0:n-1), vcat(1., zeros(n-1))) * I
R[:, 1] = b * R[:, 1]
R = R ./ sqrt(tr(R * R') / n)
mu = .5
conc_param = [sqrt(mu / k)] * (Z' * Z)^(-.5) * ones(k, 1)
s2u = 1.
s2v = 1.
rho = 1.
suv = sqrt(s2u * s2v) * rho
nuvec = 3:2:5

h = 1e-6
Zp = Z * conc_param
Pz = Z * inv(Z' * Z) * Z'
xvec = beta-3:.02:beta+3
#E = eigen(.5*Pz+.5*Pz')
#[pp,ll]=eig(.5*Pz+.5*Pz')
nx = length(xvec)
pdf = zeros(length(nuvec), nx)
cdf = similar(pdf)
cdf2 = similar(pdf)
pm = similar(pdf)
spapdf = similar(pdf)
spacdf = similar(pdf)
spacdf2 = similar(pdf)
spapm =  similar(pdf)
#Threads.@threads
for nuloop=1:length(nuvec)
    nu = nuvec[nuloop]
    for loop = 1 : nx
        Si = ((nu-2) / nu * [s2u suv; suv s2v])^.5 # rescale to make Si^2 the covariance matrix
        S = kron(Si, R*R')
        b = xvec[loop] - beta
        a0 = (-b * (Zp' * Zp))[1]
        a = [Zp; -2*b*Zp]
        A = .5*[Pz*0 Pz; Pz -2*b*Pz]
        SAS = S * A * S
        SAS = .5 * (SAS + SAS')
        E = eigen(SAS)
        P = E.vectors
        omega = E.values
        d = a' * S * P
        ccdf, p = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=false)
        cdf[nuloop, loop] = 1.0 - ccdf
        pm[nuloop, loop] = p
        ccdf, p = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=true, order=order)
        spacdf[nuloop, loop] = 1.0 - ccdf
        spapm[nuloop, loop] = p
    end
    xvech = xvec .+ h
    for loop = 1 : nx
        Si = ((nu-2) / nu * [s2u suv; suv s2v])^.5 # rescale to make Si^2 the covariance matrix
        S = kron(Si, R*R')
        b = xvech[loop] - beta
        a0 = (-b * (Zp' * Zp))[1]
        a = [Zp; -2*b*Zp]
        A = .5*[Pz*0 Pz; Pz -2*b*Pz]
        SAS = S * A * S
        SAS = .5*(SAS+SAS')
        E = eigen(SAS)
        P = E.vectors
        omega = E.values
        d = a' * S * P
        ccdf, p = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=false)
        cdf2[nuloop, loop] = 1.0 - ccdf
        pm[nuloop, loop] = p
        ccdf, p = qfmgh(0., a0, a[:], A, S, zeros(2n), zeros(2n), -nu/2, nu, 0., do_spa=true, order=order)
        spacdf2[nuloop, loop] = 1.0 - ccdf
        spapm[nuloop, loop] = p
    end
end
pdf .= (cdf2 .- cdf) ./ h
spapdf .= (spacdf2 .- spacdf) ./ h

spapdf[abs.(pdf.-spapdf) .> .1] .= NaN # filter errors due to nonsingularity around the mean
if doplot
    plot(xvec, pdf', color=collect(1:7)', labels=permutedims("v=".*string.(nuvec)))
    plot!(xvec, spapdf', color=collect(1:7)', labels="", ls=:dash)
    xlabel!("x")
    ylabel!("pdf")
    title!("Exact density (solid) and SPA (dashes)")
    savefig("pdfs.svg")
end

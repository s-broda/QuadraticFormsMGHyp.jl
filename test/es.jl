
portfolios = 1:16 # these are the portfolios from Broda (2012)
portfolios = portfolios[1]
df = 5 # degrees of freedom
lam = -df / 2
chi = float(df)
psi = 0.
m = 10 # number of factors=number of options
S0 = ones(m) * 100.0 # initial prices
K = S0 # strikes
r = .05 # risk-free rate
sigma = ones(m) * .3 # annual return volatility
put=[false, false, false, false, false, true, true, true, true, true] # true for put, false for call
h = 1 / 252 # forecast horizon of x days = x/250 years

xvec = 3.5:.1:15
ccdf = zeros(length(xvec), length(portfolios))
pm = similar(ccdf)
ccdfspa = similar(ccdf)
pmspa = similar(ccdf)
for pfloop in portfolios
    if pfloop < 9
        T = ones(m) * 126 / 252 # time to maturity in years
    else
        T = ones(m) * 21/252 # time to maturity in years
    end
    if mod(pfloop, 2) == 1
        number = -ones(m) # number of options (negative for short positions)
    else
        number = ones(m) # number of options (negative for short positions)
    end
    if pfloop == 3 || pfloop == 4 || pfloop == 7 || pfloop == 8 || pfloop == 11 || pfloop == 12 || pfloop == 15 || pfloop == 16
        deltahedged = true # include stocks to delta-hedge
    else
        deltahedged = false
    end
    if (4 < pfloop && pfloop < 9) || pfloop > 12
        equicorrelated = true # equicorrelated or uncorrelated assets
    else
        equicorrelated = false
    end

    if equicorrelated
        R = .5 * ones(m, m) + .5 * I # the correlation matrix
    else
        R = I
    end
    # get initial option prices
    price = blsprice.(S0, K, r, T, sigma, 0.0, put)
    delta = blsdelta.(S0, K, r, T, sigma, 0.0, put)
    gamma = blsgamma.(S0, K, r, T, sigma, 0.0, put)
    theta = blstheta.(S0, K, r, T, sigma, 0.0, put)
    if deltahedged
        stocks = -number .* delta
    else
        stocks = zeros(m)
    end
    a = - delta .* number - stocks
    A = diagm(0=>-.5 * gamma .* number)
    portfolio_value = sum(price .* number) + sum(stocks .* S0)
    sigs = diagm(0=>sigma .* S0 * sqrt(h)) # price volatilities 250*h days ahead
    Sigma = sigs * R * sigs * (df-2) / df # the dispersion matrix
    C = ((Sigma + Sigma') / 2) ^ .5 # and its square root
    a0 = -h * sum(number.*theta)
    c, p = qfmgh(xvec, a0, a, A, C, zeros(m), zeros(m), lam, chi, psi; do_spa=false)
    ccdf[:, pfloop] = c
    pm[:, pfloop] = p
    c, p = qfmgh(xvec, a0, a, A, C, zeros(m), zeros(m), lam, chi, psi; do_spa=true, order=2)
    ccdfspa[:, pfloop] = c
    pmspa[:, pfloop] = p
end
if doplot
    plot(ccdf[:, 1], pm[:, 1], labels="")
    plot!(ccdfspa[:, 1], pmspa[:, 1], ls=:dash, labels="")
    xlabel!("VaR level")
    ylabel!("ES")
    title!("Exact expected shortfall (solid) and SPA (dashes)")
    savefig("ES.svg")
end

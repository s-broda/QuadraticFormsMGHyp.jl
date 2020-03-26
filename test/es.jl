import QuadraticFormsMGHyp: lklam

portfolio = 1 # between 1 and 16; these are the portfolios from Broda (2012)
xvec = 3.5:.01:17.5 # results in VaR levels up to 10% for portfolio 1 for NIG with chi=psi=1

lam = -.5
chi = 1
psi = 1
m = 10 # number of factors=number of options
S0 = ones(m) * 100.0 # initial prices
K = S0 # strikes
r = .05 # risk-free rate
sigma = ones(m) * .3 # annual return volatility
put=[false, false, false, false, false, true, true, true, true, true] # true for put, false for call
h = 1 / 252 # forecast horizon of x days = x/250 years

if portfolio < 9
    T = ones(m) * 126 / 252 # time to maturity in years
else
    T = ones(m) * 21/252 # time to maturity in years
end
if mod(portfolio, 2) == 1
    number = -ones(m) # number of options (negative for short positions)
else
    number = ones(m) # number of options (negative for short positions)
end
if portfolio == 3 || portfolio == 4 || portfolio == 7 || portfolio == 8 || portfolio == 11 || portfolio == 12 || portfolio == 15 || portfolio == 16
    deltahedged = true # include stocks to delta-hedge
else
    deltahedged = false
end
if (4 < portfolio && portfolio < 9) || portfolio > 12
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
V = exp(lklam(lam+1, chi, psi) - lklam(lam, chi, psi)) # variance, assuming gamma=0
Sigma = sigs * R * sigs / V # the dispersion matrix
C = ((Sigma + Sigma') / 2) ^ .5 # and its square root
a0 = -h * sum(number.*theta)
ccdf, pm = qfmgh(xvec, a0, a, A, C, zeros(m), zeros(m), lam, chi, psi; do_spa=false)
ccdfspa, pmspa = qfmgh(xvec, a0, a, A, C, zeros(m), zeros(m), lam, chi, psi; do_spa=true, order=2)
if doplot
    plot(ccdf, pm, labels="", lw=1.25)
    plot!(ccdfspa, pmspa, ls=:dash, lw=1.25, labels="")
    xlabel!("VaR level \$\\alpha\$")
    ylabel!("\$ES_L^{(\\alpha)}\$")
    title!("Exact (solid) and approximate (dashes) expected shortfall")
    savefig("ES.pdf")
    savefig("ES.eps")
end

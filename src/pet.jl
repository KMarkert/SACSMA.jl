
function hargreaves(
    forcings::DataFrame;
    tminCol::Symbol=:tmin,
    tmaxCol::Symbol=:tmax,
    dateCol::Symbol=:datetime
)

    dts = forcings[!,dateCol]
    tmin = forcings[!,tminCol]
    tmax = forcings[!,tmaxCol]
    len = length(tmax)
    Gsc = 367
    lhov = 2.257

    doy = map(dayofyear, dts)

    tavg = map(mean,zip(tmin,tmax))

    eto = zeros(len)

    for (i,t) in enumerate(doy)
        b = 2 * pi * (Int16(t)/365)
        Rav = 1.00011 + 0.034221*cos(b) + 0.00128*sin(b) + 0.000719*cos(2*b) + 0.000077*sin(2*b)
        Ho = ((Gsc * Rav) * 86400)/1e6

        eto[i] = (0.0023 * Ho * (tmax[i]-tmin[i])^0.5 * (tavg[i]+17.8))
    end

    return eto
end

function hamon(forcings::DataFrame, par::Float64, lat::Float64; tavgCol::Symbol=:tavg, dateCol::Symbol=:dataetime)

    dts = forcings[!,dateCol]
    tavg = forcings[!,tavgCol]
    len = length(tavg)

    doy = map(dayofyear, dts)

    θ = @. 0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (doy - 186)))
	pi_v = @. asin(0.39795 * cos(θ))
    daylighthr = @. 24 - 24/pi * acos((sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(pi_v))/(cos(lat * pi/180) * cos(pi_v)))

    esat = @. 0.611 * exp(17.27 * tavg/(237.3 + tavg))

    eto = @. par * 29.8 * daylighthr * (esat/(tavg + 273.2))

  return eto
end


function snow17(forcings::DataFrame, elev::Float64, par::AbstractArray; prcpCol::Symbol=:precip, tavgCol::Symbol=:tavg, dateCol::Symbol=:datetime, initstate::AbstractArray=[0.,0.,0.,0.])
    # Define constants
    dtt = 24  # time interval of temperature data
    dtp = 24  # time interval of prcipitation data

    # SET PARAMETERS
    SCF    = par[1]
    PXTEMP = par[2]
    MFMAX  = par[3]
    MFMIN  = par[4]
    UADJ   = par[5]
    MBASE  = par[6]
    TIPM   = par[7]
    PLWHC  = par[8]
    NMF    = par[9]
    DAYGM  = par[10]

    dts = forcings[!,dateCol]
    tavg = forcings[!,tavgCol]
    prcp = forcings[!,prcpCol]
    len = length(tavg)

    doy = map(dayofyear, dts)

    
    meltNrain = Array{Float64}(undef,len)

    # Set initial states
    W_i, ATI, W_q, Deficit = initstate

    # LOOP THROUGH EACH PERIOD
    for i in 1:len

        # Set current temperature and precipitation
        Ta = tavg[i]  # Air temperature at this time step [deg C]
        Pr = prcp[i]  # Precipitation at this time step [mm]

        # FORM OF PRECIPITATION
        if Ta <= PXTEMP
            # Air temperature is cold enough for snow to occur
            SNOW = Pr
            RAIN = 0.0
        else
            # Air temperature is warm enough for rain
            SNOW = 0.0
            RAIN =  Pr
        end

        # ACCUMULATION OF THE SNOW COVER
        Pn = SNOW * SCF  # Water equivalent of new snowfall [mm]
        W_i = W_i + Pn  # Water equivalent of the ice portion of the snow cover [mm]
        E = 0  # Excess liquid water in the snow cover

        # ENERGY EXCHANGE AT SNOW/AIR SURFACE DURING NON-MELT PERIODS

        # Seasonal variation in the non-rain melt factor (Assume a year has 365 days)
        N_Mar21 = doy[i] - 80

        Sv = (0.5 * sin((N_Mar21 * 2 * pi)/365)) + 0.5  # Seasonal variation
        Av = 1  # Seasonal variation adjustment, Av<-1.0 when lat < 54N
        Mf = dtt/6 * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN)  # Seasonally varying non-rain melt factor

        # New snow temperature and heat deficit from new snow
        if (Ta < 0)
            T_snow_new = Ta
        else
            T_snow_new = 0
        end

        # Change in the heat deficit due to new snowfall [mm], 80 cal/g: latent heat of fusion, 0.5 cal/g/C:
        # specific heat of ice
        delta_HD_snow = -(T_snow_new * Pn)/(80/0.5)

        # Heat Exchange due to a temperature gradient change in heat deficit due to a temperature gradient
        # [mm]
        delta_HD_T = NMF * dtp/6 * Mf/MFMAX * (ATI - T_snow_new)

        # Update ATI[Antecedent Temperature Index]
        if (Pn > 1.5 * dtp)
            ATI = T_snow_new  #Antecedent temperature index
        else
            TIPM_dtt = 1 - ((1 - TIPM)^(dtt/6))
            ATI = ATI + TIPM_dtt * (Ta - ATI)
        end

        ATI = min(ATI, 0)

        # SNOW MELT
        T_rain = max(Ta, 0)  # Temperature of rain (deg C), Ta or 0C, whichever greater

        if (RAIN > 0.25 * dtp)
            # Rain-on-Snow Melt
            stefan = 6.12 * (10^(-10))  # Stefan-Boltzman constant (mm/K/hr)
            e_sat = 2.7489 * (10^8) * exp((-4278.63/(Ta + 242.792)))  # Saturated vapor pressure at Ta (mb)
            P_atm = 33.86 * (29.9 - (0.335 * (elev/100)) + (0.00022 * ((elev/100)^2.4)))  # Atmospheric pressure (mb) where elevation is in HUNDREDS of meters (this is incorrectly stated in the manual)
            term1 = stefan * dtp * (((Ta + 273)^4) - (273^4))
            term2 = 0.0125 * RAIN * T_rain
            term3 = 8.5 * UADJ * (dtp/6) * ((0.9 * e_sat - 6.11) + (0.00057 * P_atm * Ta))
            Melt = term1 + term2 + term3
            Melt = max(Melt, 0)

        elseif ((RAIN <= 0.25 * dtp) && (Ta > MBASE))
            # Non-Rain Melt
            Melt = (Mf * (Ta - MBASE) * (dtp/dtt)) + (0.0125 * RAIN * T_rain)
            Melt = max(Melt, 0)

        else
            Melt = 0
        end

        # Ripeness of the snow cover W_i : water equivalent of the ice portion of the snow cover W_q :
        # liquide water held by the snow W_qx: liquid water storage capacity Qw : Amount of available water
        # due to melt and rain

        Deficit = max(Deficit + delta_HD_snow + delta_HD_T, 0)  # Deficit <- heat deficit [mm]

        if (Deficit > (0.33 * W_i))
            # limits of heat deficit
            Deficit = 0.33 * W_i
        end

        if (Melt < W_i)
            W_i = W_i - Melt
            Qw = Melt + RAIN
            W_qx = PLWHC * W_i

            if ((Qw + W_q) > (Deficit + Deficit * PLWHC + W_qx))
                # THEN the snow is RIPE

                E = Qw + W_q - W_qx - Deficit - (Deficit * PLWHC)  # Excess liquid water [mm]
                W_i = W_i + Deficit  # W_i increases because water refreezes as heat deficit is decreased
                W_q = W_qx + PLWHC * Deficit  # fills liquid water capacity
                Deficit = 0

            elseif ((Qw + W_q) >= Deficit)

                # & [[Qw + W_q] <= [[Deficit*[1+PLWHC]] + W_qx]] THEN the snow is NOT yet ripe, but ice is being
                # melted

                E = 0
                W_i = W_i + Deficit  # W_i increases because water refreezes as heat deficit is decreased
                W_q = W_q + Qw - Deficit
                Deficit = 0

            elseif ((Qw + W_q) < Deficit)
                # elseif [[Qw + W_q] < Deficit]

                # THEN the snow is NOT yet ripe
                E = 0
                W_i = W_i + Qw + W_q  # W_i increases because water refreezes as heat deficit is decreased
                Deficit = Deficit - Qw - W_q
            end

        else
            Melt = W_i + W_q  # Melt >= W_i
            W_i = 0
            W_q = 0
            Qw = Melt + RAIN
            E = Qw
            # SWE = 0

        end

        if (Deficit == 0)
            ATI = 0
        end

        # Constant daily amount of melt which takes place at the snow-soil interface whenever there is a
        # snow cover
        if (W_i > DAYGM)

            gmwlos = (DAYGM/W_i) * W_q
            gmslos = DAYGM
            gmro = gmwlos + gmslos
            W_i = W_i - gmslos
            W_q = W_q - gmwlos

            E = E + gmro
            SWE = W_i + W_q

        else
          gmro = W_i + W_q
          W_i = 0
          W_q = 0
          E = E + gmro
          SWE = 0
        end

        meltNrain[i] = E
    end

    return meltNrain
end

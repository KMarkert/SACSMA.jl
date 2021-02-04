
function sacsma(forcings::DataFrame, par::AbstractArray; prcpCol::Symbol=:precip, petCol::Symbol=:pet, initstate::AbstractArray=[0,0,500,500,500,0])

    uztwm  =  par[1]    # Upper zone tension water capacity [mm]
    uzfwm  =  par[2]    # Upper zone free water capacity [mm]
    lztwm  =  par[3]    # Lower zone tension water capacity [mm]
    lzfpm  =  par[4]    # Lower zone primary free water capacity [mm]
    lzfsm  =  par[5]    # Lower zone supplementary free water capacity [mm]
    uzk    =  par[6]    # Upper zone free water lateral depletion rate [1/day]
    lzpk   =  par[7]    # Lower zone primary free water depletion rate [1/day]
    lzsk   =  par[8]    # Lower zone supplementary free water depletion rate [1/day]
    zperc  =  par[9]    # Percolation demand scale parameter [-]
    rexp   =  par[10]   # Percolation demand shape parameter [-]
    pfree  =  par[11]   # Percolating water split parameter (decimal fraction)
    pctim  =  par[12]   # Impervious fraction of the watershed area (decimal fraction)
    adimp  =  par[13]   # Additional impervious areas (decimal fraction)
    riva   =  par[14]   # Riparian vegetation area (decimal fraction)
    side   =  par[15]   # The ratio of deep recharge to channel base flow [-]
    rserv  =  par[16]   # Fraction of lower zone free water not transferrable (decimal fraction)

    # Initial Storage States (SAC-SMA)
    uztwc = initstate[1] # Upper zone tension water storage
    uzfwc = initstate[2] # Upper zone free water storage
    lztwc = initstate[3] # Lower zone tension water storage
    lzfsc = initstate[4] # Lower zone supplementary free water storage
    lzfpc = initstate[5] # Upper zone primary free water storage
    adimc = initstate[6] # Additional impervious area storage

    pet = forcings[!,petCol]
    prcp = forcings[!,prcpCol]
    len = length(pet)

    # RESERVOIR STATE ARRAY INITIALIZATION
    simflow   = Array{Float64}(undef,len) # total outflow
    base_tot  = Array{Float64}(undef,len) # base flow
    surf_tot  = Array{Float64}(undef,len) # surface runoff

    thres_zero  = 0.00001 # Threshold to be considered as zero
    parea       = 1 - adimp - pctim

    for i in 1:len

        ### Set input precipitation and potential evapotranspiration
        pr = prcp[i]
        edmnd = pet[i]

        ## Compute for different compnents...
        # ET(1), ET from Upper zone tension water storage
        et1 = edmnd * uztwc/uztwm
        red = edmnd - et1  # residual ET demand
        uztwc = uztwc - et1

        # ET(2), ET from upper zone free water storage
        et2 = 0

        # in case et1 > uztws, no water in the upper tension water storage
        if (uztwc <= 0)
            et1 = et1 + uztwc #et1 = uztwc
            uztwc = 0
            red = edmnd - et1

            # when upper zone free water content is less than residual ET
            if (uzfwc < red)

                # all content at upper zone free water zone will be gone as ET
                et2 = uzfwc
                uzfwc = 0
                red = red - et2
                uztwc = (uztwc < thres_zero) ? 0 : uztwc
                uzfwc = (uzfwc < thres_zero) ? 0 : uzfwc

            # when upper zone free water content is more than residual ET
            else
                et2 = red  # all residual ET will be gone as ET
                uzfwc = uzfwc - et2
                red = 0
            end

        # in case et1 <= uztws, all maximum et (et1) are consumed at uztwc,
        # so no et from uzfwc (et2=0)
        else

            # There's possibility that upper zone free water ratio exceeds
            #upper zone tension water ratio. If so, free water is transferred to
            #tension water storage

            if ((uztwc / uztwm) < (uzfwc / uzfwm))
                uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
                uztwc = uztwm * uzrat
                uzfwc = uzfwm * uzrat
            end

            uztwc = (uztwc < thres_zero) ? 0 : uztwc
            uzfwc = (uzfwc < thres_zero) ? 0 : uzfwc

        end

        # ET(3), ET from Lower zone tension water storage when residual ET > 0
        et3 = red * lztwc / (uztwm + lztwm) #residual ET is always bigger than ET(3)
        lztwc = lztwc - et3

        # if lztwc is less than zero, et3 cannot exceed lztws
        if(lztwc < 0)
            et3   = et3 + lztwc  # et3 = lztwc
            lztwc = 0
        end

        # Water resupply from Lower free water storages to Lower tension water storage
        saved  = rserv * (lzfpm + lzfsm)
        ratlzt = lztwc / lztwm
        ratlz  = (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved)

        # water is first taken from supplementary water storage for resupply
        if (ratlzt < ratlz)

            del = (ratlz - ratlzt) * lztwm
            lztwc = lztwc + del  # Transfer water from lzfss to lztws
            lzfsc = lzfsc - del

            # if tranfer exceeds lzfsc then remainder comes from lzfps
            if (lzfsc < 0)
                lzfpc = lzfpc + lzfsc
                lzfsc = 0
            end
        end

        lztwc = (lztwc < thres_zero) ? 0 : lztwc

        # ET(5), ET from additional impervious (ADIMP) area
        # ????? no idea where this come from, I think there's a possibility that et5 can be negative values
        et5   = et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm)
        adimc = adimc - et5

        if (adimc < 0)
            #et5 cannot exceed adims
            et5 = et5 + adimc # et5 = adimc
            adimc = 0
        end

        et5 = et5 * adimp

        # Time interval available moisture in excess of uztw requirements
        twx = pr + uztwc - uztwm

        # all moisture held in uztw- no excess
        if (twx < 0)
            uztwc = uztwc + pr
            twx = 0
            # moisture available in excess of uztw storage
        else
            uztwc = uztwm
        end

        # for now twx is excess rainfall after filling the uztwc
        adimc = adimc + pr - twx

        # Compute Impervious Area Runoff
        roimp = pr * pctim

        # Initialize time interval sums
        sbf   = 0  # Sum of total baseflow(from primary and supplemental storages)
        ssur  = 0  # Sum of surface runoff
        sif   = 0  # Sum of interflow
        sperc = 0  # Time interval summation of percolation
        sdro  = 0  # Sum of direct runoff from the additional impervious area

        # Determine computational time increments for the basic time interval
        ninc = floor(1.0 + 0.2*(uzfwc+twx))  # Number of time increments that interval is divided into for further soil-moisture accountng
        dinc = 1.0 / ninc                    # Length of each increment in days
        pinc = twx / ninc                    # Amount of available moisture for each increment

        # Compute free water depletion fractions for the time increment
        #(basic depletions are for one day)
        duz   = 1 - (1 - uzk)^dinc
        dlzp  = 1 - (1 - lzpk)^dinc
        dlzs  = 1 - (1 - lzsk)^dinc

        # Start incremental for-loop for the time interval
        for n in 1:ninc

            adsur = 0 # Amount of surface runoff. This will be updated.

            # Compute direct runoff from adimp area
            ratio = (adimc - uztwc) / lztwm
            ratio = (ratio < 0) ? 0 : ratio

            # Amount of direct runoff from the additional impervious area
            addro = pinc*(ratio^2)

            # Compute baseflow and keep track of time interval sum
            # Baseflow from free water primary storage
            bf_p = lzfpc * dlzp
            lzfpc = lzfpc - bf_p

            if (lzfpc <= 0.0001)
                bf_p  = bf_p + lzfpc
                lzfpc = 0
            end

            sbf = sbf + bf_p

            # Baseflow from free water supplemental storage
            bf_s  = lzfsc * dlzs
            lzfsc = lzfsc - bf_s

            if (lzfsc <= 0.0001)
                bf_s = bf_s + lzfsc
                lzfsc = 0
            end

            # Total Baseflow from primary and supplemental storages
            sbf = sbf + bf_s

            # Compute PERCOLATION- if no water available then skip.
            if ((pinc + uzfwc) <= 0.01)
                uzfwc = uzfwc + pinc

            else

                # Limiting drainage rate from the combined saturated lower zone storages
                percm = lzfpm * dlzp + lzfsm * dlzs
                perc = percm * uzfwc / uzfwm

                # DEFR is the lower zone moisture deficiency ratio
                defr = 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm)

                defr = (defr < 0) ? 0 : defr

                perc = perc * (1.0 + zperc * (defr^rexp))

                # Note. . . percolation occurs from uzfws before pav is added

                # Percolation rate exceeds uzfws
                perc = (perc >= uzfwc) ? uzfwc : perc

                uzfwc = uzfwc - perc    # Percolation rate is less than uzfws.

                # Check to see if percolation exceeds lower zone deficiency.
                check = lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm
                if(check > 0)
                    perc = perc - check
                    uzfwc = uzfwc + check
                end

                # SPERC is the time interval summation of PERC
                sperc = sperc + perc

                # Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
                del = uzfwc * duz # The amount of interflow
                sif = sif + del
                uzfwc = uzfwc - del

                # Distribute percolated water into the lower zones. Tension water
                # must be filled first except for the PFREE area. PERCT is
                # percolation to tension water and PERCF is percolation going to
                # free water.

                perct = perc * (1.0 - pfree)  # Percolation going to the tension water storage
                if ((perct + lztwc) <= lztwm)

                    lztwc = lztwc + perct
                    percf = 0 # Pecolation going to th lower zone free water storages

                else
                    percf = lztwc + perct - lztwm
                    lztwc = lztwm
                end

                # Distribute percolation in excess of tension requirements among the free water storages.
                percf = percf + (perc * pfree)
                if(percf != 0)

                    # Relative size of the primary storage as compared with total lower zone free water storages.
                    hpl = lzfpm / (lzfpm + lzfsm)

                    # Relative fullness of each storage.
                    ratlp = lzfpc / lzfpm
                    ratls = lzfsc / lzfsm

                    # The fraction going to primary
                    fracp = hpl * 2 * (1 - ratlp) / (2 - ratlp - ratls)

                    fracp = (fracp > 1.0) ? 1.0 : fracp

                    percp = percf * fracp # Amount of the excess percolation going to primary
                    percs = percf - percp # Amount of the excess percolation going to supplemental
                    lzfsc = lzfsc + percs


                    if (lzfsc > lzfsm)
                        percs = percs - lzfsc + lzfsm
                        lzfsc = lzfsm
                    end

                    lzfpc = lzfpc + percf - percs

                    # Check to make sure lzfps does not exceed lzfpm
                    if (lzfpc >= lzfpm)
                        excess = lzfpc - lzfpm
                        lztwc = lztwc + excess
                        lzfpc = lzfpm
                    end
                end


                # Distribute PINC between uzfws and surface runoff
                if (pinc != 0)

                    # check if pinc exceeds uzfwm
                    if((pinc + uzfwc) <= uzfwm)

                        uzfwc = uzfwc + pinc  # no surface runoff
                    else
                        sur = pinc + uzfwc - uzfwm # Surface runoff
                        uzfwc = uzfwm

                        ssur = ssur + (sur * parea)

                        # ADSUR is the amount of surface runoff which comes from
                        # that portion of adimp which is not currently generating
                        # direct runoff. ADDRO/PINC is the fraction of adimp
                        # currently generating direct runoff.
                        adsur = sur * (1.0 - addro / pinc)
                        ssur = ssur + adsur * adimp

                    end
                end
            end

            adimc = adimc + pinc - addro - adsur
            if (adimc > (uztwm + lztwm))
                addro = addro + adimc - (uztwm + lztwm)
                adimc = uztwm + lztwm
            end

            # Direct runoff from the additional impervious area
            sdro  = sdro + (addro * adimp)

            adimc = (adimc < thres_zero) ? 0 : adimc

        end # END of incremental for loop

        # Compute sums and adjust runoff amounts by the area over which they are generated.

        # EUSED is the ET from PAREA which is 1.0 - adimp - pctim
        eused = et1 + et2 + et3
        sif = sif * parea

        # Separate channel component of baseflow from the non-channel component
        tbf = sbf * parea   # TBF is the total baseflow
        bfcc = tbf / (1 + side)    # BFCC is baseflow, channel component

        # Ground flow and Surface flow
        base = bfcc                       # Baseflow and Interflow are considered as Ground inflow to the channel
        surf = roimp + sdro + ssur + sif  # Surface flow consists of Direct runoff and Surface inflow to the channel

        # ET(4)- ET from riparian vegetation.
        et4 = (edmnd - eused) * riva  # no effect if riva is set to zero

        # Compute total evapotransporation - TET
        eused = eused * parea
        tet = eused + et4 + et5

        # Check that adims >= uztws
        adimc = (adimc < uztwc) ? uztwc : adimc

        # Total inflow to channel for a timestep
        tot_outflow = surf + base - et4

        ### ------- Adjustments to prevent negative flows -------------------------#

        # If total outflow <0 surface and baseflow needs to be updated
        if (tot_outflow < 0)

            tot_outflow = 0; surf = 0; base = 0;

        else

            surf_remainder = surf - et4
            surf = max(0,surf_remainder)

            if (surf_remainder < 0) # In this case, base is reduced
                base = base + surf_remainder
                base = max(base,0)
            end
        end

        # Total inflow to channel for a timestep
        simflow[i]  = tot_outflow
        surf_tot[i] = surf
        base_tot[i] = base

    end #close time-loop

    return simflow, surf_tot, base_tot

end # close function

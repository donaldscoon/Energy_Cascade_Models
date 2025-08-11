import Init
import numpy as np
import pandas as pd
import os

def mpr(PPFD, CO2):
    # Define a list of coefficient expressions
    terms = [
        (1/PPFD)*(1/CO2), (1/PPFD), (CO2/PPFD), (CO2**2/PPFD), (CO2**3/PPFD),
        (1/CO2), 1, CO2, (CO2**2), (CO2**3),
        PPFD*(1/CO2), PPFD, PPFD*CO2, PPFD*(CO2**2), PPFD*(CO2**3),
        (PPFD**2)*(1/CO2), (PPFD**2), (PPFD**2)*CO2, (PPFD**2)*(CO2**2), (PPFD**2)*(CO2**3),
        (PPFD**3)*(1/CO2), (PPFD**3), (PPFD**3)*CO2, (PPFD**3)*(CO2**2), (PPFD**3)*(CO2**3)
    ]
    
    return terms

def t_A_coefficients():
    """ Canopy Closure t_A """
    # t_A coefficients (tac) values for lettuce originate from Ewert table 4-115 
    t_A_c = [0,1.0289*(10**4),-3.7018,0 ,3.6648*(10**-7),
         0,1.7571,0,2.3127*(10**-6),0,
         1.8760,0,0,0,0,
         0,0,0,0,0,
         0,0,0,0,0]
    return t_A_c

def calc_t_A(PPFD, CO2):
    terms = mpr(PPFD, CO2)
    t_A_c = t_A_coefficients()

    t_A = []
    
    for term, t_A_c_value in zip(terms, t_A_c):
        multiplied_value = term * t_A_c_value
        t_A.append(multiplied_value)

    t_A = sum(t_A)

    return t_A
        # print(term, t_A_c_value)

def CQY_max_coefficients():
    """ Canopy Quantum Yield Equation """
    # CQY_max Coefficients for lettuce from ewert table 4-102, 
    CQY_m_c = [0, 0, 0, 0, 0, 
               0, 4.4763*(10**-2), 5.163*(10**-5), -2.075*(10**-8), 0,
               0, -1.1701*(10**-5), 0, 0, 0,
               0, 0, -1.9731*(10**-11), 8.9265*(10**-15), 0, 
               0, 0, 0, 0, 0]

    return CQY_m_c

def calc_CQY_max(PPFD, CO2):
    terms = mpr(PPFD, CO2)
    CQY_m_C = CQY_max_coefficients()

    CQY_max = []
    
    for term, CQY_m_value in zip(terms, CQY_m_C):
        multiplied_value = term * CQY_m_value
        CQY_max.append(multiplied_value)

    CQY_max = sum(CQY_max)

    return CQY_max

def save_simulation(df_records, path, model):
    os.makedirs(path, exist_ok=True)
    df_records.to_csv(f'{path}/{model}_SIM_OUT.csv', index=False, header=True)
    print(f'Simulation data saved to {path}/{model}_SIM_OUT.csv')
    return

def CAV(input_file, TCB, t_A, output, save_output=True):
    """
    Calculates the Total Edible Biomass (TEB) in response to daily environmental data.

    This function reads the input data from a given file line by line,
    initializes model parameters, and then caclulate the crop yield for lettuce at 19.2 plants/m2.

    This code recreates the work detailed in the peer-reviewed publication 
    Cavazzoni, J. Using Explanatory Crop Models to Develop Simple Tools for Advanced Life Support System Studies. 
    Advances in Space Research 2004, 34, 1528–1538, doi:10.1016/j.asr.2003.02.073.

    Args:
        input_file (str): Path to the input file containing necessary data.
        TCB (float): Initial Total Crop Biomass (TCB) value. Default = 0
        t_A (int): the number of days after seeding at which the canopy is "closed" with no space remaining between plants
        output (str): Path to output file where the simulation data will be saved.
        save_output (bool): Whether to save simulation data to disk.

    Returns:
        DataFrame: A dataframe containing all the inputs, and outputs of the simulation
    """
    inputs, t_M, t_T = Init.read_input_file(input_file)
    BCF, XFRT, OPF, g_A, A_max, t_Q, t_E, MW_W, CQY_min, CUE_max, CUE_min, D_PG, p_W, n, MWC, MW_O2, MW_CO2, a, b, WBF, DRY_FR, NC_FR = Init.model_parameters()
    _, res, ts_to_harvest, df_records, _, TEB= Init.INTIALIZATION(t_M)
    model = 'CAV'
    TEB = TCB * XFRT

    for index, row in inputs.iterrows():
        # Access data in each row
        t = row['t']
        TEMP = row['TEMP']
        RH = row['RH']
        PPFD = row['PPFD']
        CO2 = row['CO2']
        H = row ['H'] 
        P_ATM = row['P_ATM']

        # Unparameterized t_A
        t_A = calc_t_A(PPFD, CO2)
        # Parameterized t_A
        # t_A = t_A
        CQY_max = calc_CQY_max(PPFD, CO2)

        if t < ts_to_harvest:
            if t < t_A:                  # before canopy closure
                A = A_max*(t/t_A)**n         # Ewert eq 4-14
            else:                        # after canopy closure
                A = A_max                    # Ewert eq 4-14
            if t<= t_Q:                  # before onset of senescence
                CQY = CQY_max                # ewert eq 4-15
                CUE_24 = CUE_max             # ewert eq 4-16
            else: 
                CQY = CQY_max - (CQY_max - CQY_min)*((t-t_Q)/(t_M-t_Q)) # ewert eq 4-15
                CUE_24 = CUE_max - (CUE_max - CUE_min)*((t-t_Q)/(t_M-t_Q)) #ewert eq 4-16
                print("Error: Utilizing CQY and CUE values without definitions")
                break
            DCG = 0.0036*H*CUE_24*A*CQY*PPFD # ewert eq 4-17 number is related to seconds in an hour
            DOP = OPF*DCG                    # ewert eq 4-18
            CGR = MWC*(DCG/BCF)            # ewert eq 4-19 number is molecular weight of carbon
            TCB += CGR                       # ewert eq 4-20
            if t > t_E:                      # accumilate edible biomass when organ formation begins
                TEB += XFRT*CGR              # ewert eq 4-21
            VP_SAT = 0.611 * np.exp((17.4 * TEMP) / (TEMP + 239)) # assumes leaf tempp=air temp. Saturated Vapor Pressure. ewert eq 4-23 numbers likely from Monje 1998
            VP_AIR = VP_SAT*(RH/100)               # Atmo Vapor Pressure ewewrt eq 4-23
            VPD = VP_SAT - VP_AIR            # Vapor Pressure Deficit ewert eq 4-23
            P_GROSS = A*CQY*PPFD             # Gross photosynthesis ewert eq 4-24
            P_NET = (((D_PG-H)/D_PG)+((H*CUE_24)/D_PG))*P_GROSS     # Net Photosynthesis ewert eq 4-25
            g_S = (1.717*TEMP-19.96-10.54*VPD)*(P_NET/CO2)        # stomatal conductance the numbers came from monje 1998, only for planophile canopies equation from ewert 4-27
            g_C = (g_A*g_S)/(g_A+g_S)                               # canopy conductance ewert 4-26
            DTR = 3600*H*(MW_W/p_W)*g_C*(VPD/P_ATM)

            df_ts = pd.DataFrame({
                'DAS': [t],
                'TEMP': [TEMP],
                'RH': [RH],
                'PPFD': [PPFD],
                'CO2': [CO2],
                'H': [H],
                'P_ATM': [P_ATM],
                't_A': [t_A],
                't_A_%': [min(t/t_A, 1)],
                'A': [A],
                'CQY': [CQY],
                'CUE_24': [CUE_24],
                'DCG': [DCG],
                'DOP': [DOP],
                'CGR': [CGR],
                'TCB': [TCB],
                'TEB': [TEB],   
                'VP_SAT': [VP_SAT],
                'VP_AIR': [VP_AIR],
                'VPD': [VPD],
                'P_GROSS': [P_GROSS],
                'P_NET': [P_NET],
                'g_S': [g_S],
                'g_C': [g_C],
                'DTR': [DTR],
            })
            df_records = pd.concat([df_records, df_ts], ignore_index=True)
    if save_output:
        save_simulation(df_records, output, model)    
    return df_records

def BOS(input_file, TCB, t_A, output, save_output=True):
    """
    Calculates the Total Edible Biomass (TEB) in response to daily environmental data.

    This function reads the input data from a given file line by line,
    initializes model parameters, and then caclulate the crop yield for lettuce at 19.2 plants/m2.

    This code recreates the work detailed in the peer-reviewed publication 
    Boscheri, G.; Kacira, M.; Patterson, L.; Giacomelli, G.; Sadler, P.; Furfaro, R.; Lobascio, C.; Lamantea, M.; Grizzaffi, L. 
    Modified Energy Cascade Model Adapted for a Multicrop Lunar Greenhouse Prototype. 
    Advances in Space Research 2012, 50, 941–951, doi:10.1016/j.asr.2012.05.025.

    Args:
        input_file (str): Path to the input file containing necessary data.
        TCB (float): Initial Total Crop Biomass (TCB) value. Default = 0
        t_A (int): the number of days after seeding at which the canopy is "closed" with no space remaining between plants
        output (str): Path to output file where the simulation data will be saved.
        save_output (bool): Whether to save simulation data to disk.

    Returns:
        DataFrame: A dataframe containing all the inputs, and outputs of the simulation
    """
    inputs, t_M, t_T = Init.read_input_file(input_file)
    BCF, XFRT, OPF, g_A, A_max, t_Q, t_E, MW_W, CQY_min, CUE_max, CUE_min, D_PG, p_W, n, MWC, MW_O2, MW_CO2, a, b, WBF, DRY_FR, NC_FR = Init.model_parameters()
    _, res, ts_to_harvest, df_records, _, TEB= Init.INTIALIZATION(t_M)
    model = 'BOS'

    TEB = TCB * XFRT

    for index, row in inputs.iterrows():
        # Access data in each row
        t = row['t']
        TEMP = row['TEMP']
        RH = row['RH']
        PPFD = row['PPFD']
        CO2 = row['CO2']
        H = row ['H'] 
        P_ATM = row['P_ATM']

        # Unparameterized t_A
        t_A = calc_t_A(PPFD, CO2)
        # # Parameterized t_A
        # t_A = t_A
        CQY_max = calc_CQY_max(PPFD, CO2)

        if t < ts_to_harvest:
            if t < t_A:                  # before canopy closure
                A = A_max*(t/t_A)**n         # boscheri eq 5
            else:                        # after canopy closure
                A = A_max                    # boscheri eq 5
            if t<= t_Q:                  # before onset of senescence
                CQY = CQY_max                # boscheri eq 3
                CUE_24 = CUE_max             # boscheri eq 4
            elif t_M > t: 
                """For lettuce the values of CQY_min and CUE_min 
                are n/a due to the assumption that the canopy does
                not senesce before harvest. I coded them anyways, it
                makes it complete for all the other crops too. For 
                crops other than lettuce remove the break statement."""
                CQY = CQY_max - (CQY_max - CQY_min)*((t-t_Q)/(t_M-t_Q)) # boscheri eq 3
                CUE_24 = CUE_max - (CUE_max - CUE_min)*((t-t_Q)/(t_M-t_Q)) # boscheri eq 4
                print(t, "Error: Utilizing CQY and CUE values without definitions")
                break 
            HCG = 0.0036*H*CUE_24*A*CQY*PPFD      # boscheri eq 2 but removed I and multiplied alpha by 24  
            HCGR = HCG*MWC*(BCF)**(-1)       # boscheri eq 6 Dry
            ######## SEE WBF FOR ASSUMPTION #############
            HWCGR = HCGR*(1-WBF)**(-1)       # boscheri eq 7 Wet
            HOP = HCG/CUE_24*OPF*MW_O2       # boscheri eq 8
            HOC = HCG/(1-CUE_24)/CUE_24*OPF*MW_O2*H/24 # paper includes "I" with a weird notation, but can't divide by I so I removed boscheri eq 9 
            VP_SAT = 0.611 * np.exp((17.4 * TEMP) / (TEMP + 239)) # boscheri eq 12
            VPD = VP_SAT*(1-(RH/100))             # boscheri eq 12
            P_NET = A*CQY*PPFD              # boscheri eq 13
            g_S = (1.717*TEMP-19.96-10.54*VPD)*(P_NET/CO2) # boscheri unlabeled equation
            g_C = (g_A*g_S)*(g_A+g_S)**(-1) # boscheri unlabeled equation
            # HTR = b*MW_W*g_C*(VPD/P_ATM)    # boscheir eq 10 Which is believed to be missing the density of water
            HTR = 3600*H*(MW_W/p_W)*g_C*(VPD/P_ATM) # similar to the original but takes into account the density of water conversion factor corrected for daily
            HCO2C = HOP*MW_CO2*MW_O2**(-1)  # boscheri eq 14
            HCO2P = HOC*MW_CO2*MW_O2**(-1)  # boscheri eq 15 
            HNC = HCGR*DRY_FR*NC_FR         # boscheri eq unlabeled
            HWC = HTR+HOP+HCO2P+HWCGR-HOC-HCO2C-HNC # boscheri eq 16
            TCB += HCGR                       # ewert eq 4-20
            if t > t_E:                      # accumilate edible biomass when organ formation begins
                TEB += XFRT*HCGR              # ewert eq 4-21
            df_ts = pd.DataFrame({
                'DAS': [t],
                'TEMP': [TEMP],
                'RH': [RH],
                'PPFD': [PPFD],
                'CO2': [CO2],
                'H': [H],
                'P_ATM': [P_ATM],
                't_A': [t_A],
                't_A_%': [min(t/t_A, 1)],
                'A': [A],
                'CQY': [CQY],
                'CUE_24': [CUE_24],
                'DCG': [HCG],
                'CGR': [HCGR],
                'WCGR': [HWCGR],
                'DOP': [HOP],
                'DOC': [HOC],
                'VP_SAT': [VP_SAT],
                'VPD': [VPD],
                'P_NET': [P_NET],
                'g_S': [g_S],
                'g_C': [g_C],
                'DTR': [HTR],
                'DCO2C': [HCO2C],
                'DCO2P': [HCO2P],
                'DNC': [HNC],
                'DWC': [HWC],
                'TCB': [TCB],
                'TEB': [TEB]
                })
            df_records = pd.concat([df_records, df_ts], ignore_index=True)
    if save_output:
        save_simulation(df_records, output, model)  
    return df_records

def AMI(input_file, TCB, t_A, output, save_output=True):
    """
    Calculates the Total Edible Biomass (TEB) and Transpiration Rate (DTR) in response to daily environmental data.

    This function reads the input data from a given file line by line,
    initializes model parameters, and then caclulate the crop yield for lettuce at 19.2 plants/m2.

    This code recreates the work detailed in the peer-reviewed publication 
    Amitrano, C.; Chirico, G.B.; De Pascale, S.; Rouphael, Y.; De Micco, V. 
    Crop Management in Controlled Environment Agriculture (CEA) Systems Using Predictive Mathematical Models. 
    Sensors 2020, 20, 3110, doi:10.3390/s20113110.
    
    Args:
        input_file (str): Path to the input file containing necessary data.
        TCB (float): Initial Total Crop Biomass (TCB) value. Default = 0
        t_A (int): the number of days after seeding at which the canopy is "closed" with no space remaining between plants
        output (str): Path to output file where the simulation data will be saved.
        save_output (bool): Whether to save simulation data to disk.

    Returns:
        DataFrame: A dataframe containing all the inputs, and outputs of the simulation
    """

    inputs, t_M, t_T = Init.read_input_file(input_file)
    BCF, XFRT, OPF, g_A, A_max, t_Q, t_E, MW_W, CQY_min, CUE_max, CUE_min, D_PG, p_W, n, MWC, MW_O2, MW_CO2, a, b, WBF, DRY_FR, NC_FR = Init.model_parameters()
    _, res, ts_to_harvest, df_records, _, TEB= Init.INTIALIZATION(t_M)
    amin_GN, amin_GON, amax_GN, amax_GON, bmin_GN, bmin_GON, bmax_GN, bmax_GON = Init.AMI_params()
    model = 'AMI'

    t_D = 1
    TEB = TCB * XFRT

    for index, row in inputs.iterrows():
        # Access data in each row
        t = row['t']
        TEMP = row['TEMP']
        RH = row['RH']
        PPFD = row['PPFD']
        CO2 = row['CO2']
        H = row ['H'] 
        P_ATM = row['P_ATM']

        if t < ts_to_harvest:
            VP_SAT = 0.611 * np.exp((17.4 * TEMP) / (TEMP + 239)) # Same as ewert and cavazzoni, though likely from Monje 1998
            VP_AIR = VP_SAT*(RH/100)                     # Same as ewert and cavazzoni, though likely from Monje 1998
            VPD = VP_SAT*(1-(RH/100))                    # Same as ewert and cavazzoni, though likely from Monje 1998
            if VPD <= 1.225:                         # Assuming the nominal inflection point as it is halfway between 0.69(NOM) and 1.76(OFFNOM)
                amin = amin_GN
                amax = amax_GN
                bmin = bmin_GN
                bmax = bmax_GN
            else:
                amin = amin_GON
                amax = amax_GON
                bmin = bmin_GON
                bmax = bmax_GON 
            if t<= t_D:                            # if timestep is before formation of edible organs
                ALPHA = amin                    # amitrano 2020 eq 15
                BETA = bmin                     # amitrano 2020 eq 15
            elif t <= t_M:                        # if timestep is after organ formation but before maturity
                ALPHA = amin+(amax-amin)*(t-t_D)/(t_M-t_D)        # amitrano 2020 eq 15
                BETA = bmin+(bmax-bmin)*(t-t_D)/(t_M-t_D)         # amitrano 2020 eq 15
            else:                                  # all other timesteps
                ALPHA = amax                    # amitrano 2020 eq 15
                BETA = bmax     
            DCG = 0.0036*H*ALPHA*PPFD              # amitrano 2020 eq 4
            DOP = OPF*DCG                          # amitrano 2020 eq 5
            CGR = MWC*(DCG/BCF)                      # amitrano 2020 eq 6
            if t > t_E:                            # if edible organ formation has begun
                TEB = CGR+TEB                      # Amitrano 2020 GN excel column I
                TCB = TEB*1.05
            P_GROSS = BETA*PPFD                    # amitrano 2020 eq 8
            P_NET = (H*ALPHA/24+BETA*(24-H)/24)*PPFD    # Amitrano 2020 eq 9
            g_S = ((1.717*TEMP)-19.96-(10.54*VPD))*(P_NET/CO2) # Amitrano 2020 eq 10 (with some nice parenthesis that don't change anything)
            g_C = g_A*g_S/(g_A+g_S)                # Amitrano 2020 eq 10
            DTR = 3600*H*(MW_W/p_W)*g_C*(VPD/P_ATM)
            df_ts = pd.DataFrame({
                'DAS': [t],
                'TEMP': [TEMP],
                'RH': [RH],
                'PPFD': [PPFD],
                'CO2': [CO2],
                'H': [H],
                'P_ATM': [P_ATM],
                't_A': [t_A],
                't_A_%': [min(t/t_A, 1)],
                'VP_SAT': [VP_SAT],
                'VP_AIR': [VP_AIR],
                'VPD': [VPD],
                'ALPHA': [ALPHA],
                'BETA': [BETA],
                'DCG': [DCG],
                'DOP': [DOP],
                'CGR': [CGR],
                'TEB': [TEB],
                'TCB': [TCB],
                'P_GROSS': [P_GROSS],
                'P_NET': [P_NET],
                'g_S': [g_S],
                'g_C': [g_C],
                'DTR': [DTR],
                })
            df_records = pd.concat([df_records, df_ts], ignore_index=True)
    if save_output:
        save_simulation(df_records, output, model)  
    return df_records

def EC(input_file, TCB, t_A, output, save_output=True):
    """
    Calculates the Total Edible Biomass (TEB) in response to daily environmental data.

    This function reads the input data from a given file line by line,
    initializes model parameters, and then caclulate the crop yield for lettuce at 19.2 plants/m2.

    This code recreates the work detailed in the peer-reviewed publication 
    Volk, T.; Bugbee, B.; Wheeler, R.M. 
    An Approach to Crop Modeling with the Energy Cascade. 
    Life Support Biosph Sci 1995, 1, 119–127.
    
    Args:
        input_file (str): Path to the input file containing necessary data.
        TCB (float): Initial Total Crop Biomass (TCB) value. Default = 0
        t_A (int): the number of days after seeding at which the canopy is "closed" with no space remaining between plants
        output (str): Path to output file where the simulation data will be saved.
        save_output (bool): Whether to save simulation data to disk.

    Returns:
        DataFrame: A dataframe containing all the inputs, and outputs of the simulation
    """
    inputs, t_M, t_T = Init.read_input_file(input_file)
    BCF, XFRT, OPF, g_A, A_max, t_Q, t_E, MW_W, CQY_min, CUE_max, CUE_min, D_PG, p_W, n, MWC, MW_O2, MW_CO2, a, b, WBF, DRY_FR, NC_FR = Init.model_parameters()
    _, res, ts_to_harvest, df_records, _, TEB= Init.INTIALIZATION(t_M)
    model = 'EC'

    TEB = TCB * XFRT

    H = 18 # Photoperiod
    T = 20 # Temperature listed but unused in Table Volk et al. 1995 Table 2
    K = 0.098 # Conversion constant Volk et al. 1995 table 2 Has plant values and time conversions in it.
    C = 0.625 # Carbon Use Efficency Volk et al. 1995 Matches value from bos and cav for lettuce

    # t_a = 22 # Time (days) of canopy closure Volk et al. 1995 Table 3 adjusted to be similar CAV/BOS
    t_a = t_A
    t_Q = 62 # Time (days) of onset of canopy senescence Volk et al. 1995 Table 3

    Q_min = 0.0125 # canopy quantum yield at t=t_m Volk et al. 1995 Table 3
    Q_max = 0.054 # canopy quantum yield until t=t_Q similar to results from cav and bos models

    for index, row in inputs.iterrows():
        # Access data in each row
        t = row['t']
        PPFD = row['PPFD']
        H = row ['H'] 

        if t <= t_a: # Before canopy closure
            A = (A_max/t_a)*t # Volk et al. 1995 Equation 6a
        elif t > t_a: # After canopy closure
            A = A_max # Volk et al. 1995 Equation 6b
        if t <= t_Q:
            Q = Q_max # Volk et al. 1995 Equation 7a
        elif t > t_Q:
            Q = Q_max - ((Q_max-Q_min)/(t_M-t_Q))*(t-t_Q) # Volk et al. 1995 Equation 7b

        P_GROSS = Q*A*PPFD # Gross Photosynthesis Volk et al. 1995 equation 1
        P_NET = C*Q*A*PPFD # Net Photosysnthesis Volk et al. 1995 equation 2
        R = (1-C)*Q*A*PPFD # Respiration Volk et al. 1995 equation 3
        CGR = K*(H*P_NET-(24-H)*R) # Crop Growth Rate Volk et al. 1995 equation 4
        TCB += CGR # Biomass Volk et al. 1995 equation 5

        dfts = pd.DataFrame({
            'DAS': [t],
            't_A': [t_A],
            't_A_%': [min(t/t_A, 1)],
            'A': [A],
            'CQY': [Q],
            'CUE_24': [C],
            'P_GROSS': [P_GROSS],
            'P_NET': [P_NET],
            'RESP': [R],
            'CGR': [CGR],
            'TCB': [TCB],
            'H': [H],
            'PPFD': [PPFD],
        })
        df_records = pd.concat([df_records, dfts], ignore_index=True) # this adds the timestep dataframe to the historical values dataframe
    if save_output:
        save_simulation(df_records, output, model)  
    return df_records
import pandas as pd

def read_input_file(file_path):
    """
    Reads input CSV file and validates required columns for crop simulation.

    Parameters:
        file_path (str): Path to the CSV data file.

    Returns:
        df (pd.DataFrame): DataFrame containing input data.
        t_M (int): Model duration, seeds + timesteps.
        t_T (int): Transplant day offset.

    Raises:
        ValueError: If required columns are missing or have invalid types.
        FileNotFoundError: If the file does not exist.
        pd.errors.ParserError: If the file cannot be parsed as CSV.
    """

    # Checking for required inputs
    required_columns = ['t', 'TEMP', 'RH', 'PPFD', 'CO2', 'H', 'P_ATM']
    try:
        df = pd.read_csv(file_path, header=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file not found: {file_path}")
    except pd.errors.ParserError as e:
        raise ValueError(f"Could not parse CSV file: {e}")

    # Returns any missing columns
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Basic type checking, e.g. all numeric except 't'
    for col in required_columns:
        if col != 't' and not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(f"Column {col} must be numeric")
        
    t_M = len(df) + 14 #transplants at 14 days after seeding
    t_T = t_M - len(df)

    return df, t_M, t_T

def INTIALIZATION(t_M):
    t = 0               # time in days
    res = 1             # model resolution (in days)
    ts_to_harvest = int(t_M/res)            # calcs the timesteps needed to set up the matrix for each ts
    df_records = pd.DataFrame({})           # creates a dataframe to store the results
    TCB = 0                                 # starting crop biomass
    TEB = 0                                 # starting total edible biomass

    return t, res, ts_to_harvest, df_records, TCB, TEB

def model_parameters():
    '''These are only the lettuce parameters'''
    BCF = 0.40          # Biomass carbon fraction ewert table 4-113
    XFRT = 0.95         # edible biomass fraction ewert table 4-112
    OPF = 1.08          # Oxygen production fraction ewert table 4-113
    g_A = 2.5           # atmospheric aerodynamic conductance ewert eq 4-27 no citations
    A_max = 0.93        # maximum fraction of PPF Absorbtion ewert pg 180
    t_Q = 50            # onset of senescence placeholder value ewert table 4-112
    t_E = 1             # time at onset of organ formation ewert table 4-112
    MW_W = 18.015       # Molecular weight of water, ewert table 4-110
    MWC = 12.0107       # molecular weight of carbon boscheri table 4
    MW_O2 = 31.9988     # molecular weight of O2 boscheri table 4
    MW_CO2 = 44.010     # molecular weight of CO2 boscheri table 4
    CQY_min = 0         # N/A minimum canopy quantum yield ewert table 4-99
    CUE_max = 0.625     # maximum carbon use efficiency ewert table 4-99
    CUE_min = 0         # N/A minimum carbon use efficiency ewert table 4-99
    D_PG = 24           # the plants diurnal cycle length assumed 24 in cavazzoni 2001
    p_W = 998.23        # density of water at 20 C, ewert table 4-110
    n = 2.5             # Ewert table 4-97 crop specific
    a = 0.0036          # boscheri table 4 similar to others but in 'a'
    b = 3600            # boscheri table 4
    WBF = XFRT          # Boscheri doesn't define this, I'm assuming that its the same as XFRT
    DRY_FR = 6.57/131.35# dry over wet biomass fraction Hanford 2004 Table 4.2.7, with part from wheeler 2003
    NC_FR = 0.034       # Hanford 2004 table 4.2.10, ugh boscheri just state the number

    return BCF, XFRT, OPF, g_A, A_max, t_Q, t_E, MW_W, CQY_min, CUE_max, CUE_min, D_PG, p_W, n, MWC, MW_O2, MW_CO2, a, b, WBF, DRY_FR, NC_FR

def AMI_params():
    '''Only the parameters for green lettuce'''
    amin_GN = 0.00691867456539118    # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    amin_GON = 0.00342717997911672   # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    amax_GN = 0.017148682744336      # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    amax_GON = 0.00952341360955465   # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    bmin_GN = 0                      # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    bmin_GON = 0.0486455477321762    # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    bmax_GN = 0.0451765692503675     # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    bmax_GON = 0.0564626043274799    # amitrano 2020 calibrated with growth chamber experiment exact value from excel
    
    return amin_GN, amin_GON, amax_GN, amax_GON, bmin_GN, bmin_GON, bmax_GN, bmax_GON

def seed():
    # This function is used to input generic "nominal" conditions for consistent results
    # Use while developing and testing new model components
    TEMP = 22
    RH = 85
    PPFD = 350
    CO2 = 500
    H = 12
    t_M = 30
    P_ATM = 101

    return TEMP, RH, PPFD, CO2, H, t_M, P_ATM

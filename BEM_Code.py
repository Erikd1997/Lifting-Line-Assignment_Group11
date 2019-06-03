import numpy as np

def CTfunction(a, glauert = False):
    """
    This function calculates the thrust coefficient as a function of induction factor 'a'
    'glauert' defines if the Glauert correction for heavily loaded rotors should be used; default value is false
    """
    CT = np.zeros(np.shape(a))
    CT = 4*a*(1-a)
    if glauert:
        CT1 = 1.1816;
        a1 = 1-np.sqrt(CT1)/2;
        CT[a>=a1] = CT1-4*(np.sqrt(CT1)-1)*(1-a[a>=a1])
    return CT
    
def ainduction(CT, glauert):
    """
    This function calculates the induction factor 'a' as a function of thrust coefficient CT 
    including Glauert's correction
    """
    a = np.zeros(np.shape(CT))
    CT1=1.816;
    CT2=2*np.sqrt(CT1)-CT1
    if glauert:
        a[CT>=CT2] = 1 + (CT[CT>=CT2]-CT1)/(4*(np.sqrt(CT1)-1))
        a[CT<CT2] = 0.5-0.5*np.sqrt(1-CT[CT<CT2])
    elif not glauert:
        a = 0.5-0.5*np.sqrt(1-CT)
    return a

def PrandtlTipRootCorrection(flowtype, r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    """
    This function calculate step combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
    given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number of blades NBlades and the axial induction factor
    """
    if flowtype == 'turbine':
        a = (1-axial_induction)
    elif flowtype == 'propeller':
        a = (1-axial_induction)
        
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*np.sqrt(1+((TSR*r_R)**2)/(a**2))
    Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0
    
    #print("Ftip",Ftip)
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*np.sqrt(1+((TSR*r_R)**2)/(a**2))
    Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    return Froot*Ftip, Ftip, Froot

def OptimalAlpha(polar_alpha, polar_cl, polar_cd):
    alpha_range = np.arange(-10,25,0.1)
    L_D = np.zeros((len(alpha_range),1))
    for i in range(len(alpha_range)):
        alpha = alpha_range[i]
        
        cl = np.interp(alpha, polar_alpha, polar_cl)
        cd = np.interp(alpha, polar_alpha, polar_cd)
        
        L_D[i] = cl/cd
    L_D_lst = L_D.tolist()
    index = L_D_lst.index(max(L_D_lst))
    Alpha_opt = alpha_range[index]
    return Alpha_opt

def OptimalChord(inflowangle, fnorm, vnorm, vtan, rho, twist, polar_alpha, polar_cl, polar_cd):
    vmag2 = vnorm**2 + vtan**2
    alpha = inflowangle*180/np.pi+twist
    a = np.cos(inflowangle)
    b = np.sin(inflowangle)
    
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    
    lift_chord = 0.5*cl*rho*vmag2
    drag_chord = 0.5*cd*rho*vmag2
    Chord_opt = fnorm/(lift_chord*a+drag_chord*b)

    return Chord_opt

# define function to determine load in the blade element
def loadBladeElement(rho, vnorm, vtan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd):
    """
    calculates the load in the blade element
    """
    vmag2 = vnorm**2 + vtan**2
    inflowangle = abs(np.arctan2(vnorm,vtan))
    alpha = inflowangle*180/np.pi+twist
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    lift = 0.5*rho*vmag2*cl*chord
    drag = 0.5*rho*vmag2*cd*chord
    fnorm = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
    ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
    gamma = 0.5*np.sqrt(vmag2)*cl*chord
    ct = ftan/(0.5*rho*vmag2*chord)
    cn = fnorm/(0.5*rho*vmag2*chord)
    return fnorm, ftan, gamma, alpha, inflowangle, ct, cn

def solveStreamtube(flowtype, rho, Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, tsr):
    """
    solve balance of momentum between blade element load and loading in the streamtube
    input variables:
    Uinf - wind speed at infinity
    r1_R,r2_R - edges of blade element, in fraction of Radius ;
    rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
    Radius is the rotor radius
    Omega -rotational velocity
    NBlades - number of blades in rotor
    """
    Area = np.pi*((r2_R*Radius)**2-(r1_R*Radius)**2) #  area streamtube
    r_R = (r1_R+r2_R)/2 # centroid
    
    # initialize variables
    a = 0.0 # axial induction
    aline = 0.0 # tangential induction factor
    
    Niterations = 100
    Erroriterations = 0.00001 # error limit for iteration process, in absolute value of induction
    
    for i in range(Niterations):
        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate velocity and loads at blade element"
        # ///////////////////////////////////////////////////////////////////////
        if flowtype == 'turbine':
            Urotor = Uinf*(1-a) # axial velocity at rotor (turbine case)
            Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
        elif flowtype == 'propeller':
            Urotor = Uinf*(1+a) # axial velocity at rotor (propeller case)
            Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
        
        
        # calculate loads in blade segment in 2D (N/m)
        fnorm, ftan, gamma, alpha, inflowangle, ct, cn = loadBladeElement(rho, Urotor, Utan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd)
        load3Daxial = fnorm*Radius*(r2_R-r1_R)*NBlades # 3D force in axial direction
        
        # ///////////////////////////////////////////////////////////////////////
        # //the block "Calculate velocity and loads at blade element" is done
        # ///////////////////////////////////////////////////////////////////////

        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////
        # // calculate thrust coefficient at the streamtube 
        CT = load3Daxial/(0.5*Area*rho*Uinf**2)
        CQ = 4*aline*(1-a)*tsr*r_R
        
        # calculate new axial induction, accounting for Glauert's correction
        if flowtype == 'turbine':
            anew = ainduction(CT, True)
        elif flowtype == 'propeller':
            anew = ainduction(CT, True)
        
        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(flowtype, r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
        #print(Prandtl)
        if (Prandtl < 0.0001):
            
            Prandtl = 0.0001 # avoid divide by zero
            #print(Prandtl)
        anew = anew/Prandtl # correct estimate of axial induction
        # calculate azimuthal induction
        aline = ftan*NBlades/(rho*2*np.pi*Urotor*Omega*2*(r_R*Radius)**2)
        aline = aline/Prandtl # correct estimate of azimuthal induction with Prandtl's correction
        # ///////////////////////////////////////////////////////////////////////////
        # // end of the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////
        
        #// test convergence of solution, by checking convergence of axial induction
        if (np.abs(a-anew) < Erroriterations): 
#            print("iterations")
#            print(i)
            break
        
        a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction
    return [a , aline, r_R, fnorm, ftan, gamma, alpha, inflowangle, ct, cn, Prandtl, Prandtltip, Prandtlroot, CQ, Urotor, Utan]


def Optimise(flowtype, fnormi, Uinf, r1_R, r2_R, Radius, rootradius_R, tipradius_R, Omega, NBlades, rho, polar_alpha, polar_cl, polar_cd, tsr):
    r_R = (r1_R+r2_R)/2 # centroid
    
    a = 1/3
    aline = a*(1-a)/(tsr**2*r_R**2)
    inflowangle = np.arctan2(2/3,tsr*r_R*(1+2/(3*tsr**2*r_R**2)))
    alpha = OptimalAlpha(polar_alpha, polar_cl, polar_cd)
    twist = alpha-inflowangle*180/np.pi
    
    Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(flowtype, r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, a);
    a = a/Prandtl
    aline = aline/Prandtl
    if flowtype == 'turbine':
        Urotor = Uinf*(1-a)# axial velocity at rotor (turbine case)
        Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
    elif flowtype == 'propeller':
        Urotor = Uinf*(1+a) # axial velocity at rotor (propeller case)
        Utan = (1-aline)*Omega*r_R*Radius # tangential velocity at rotor
    
    chord = OptimalChord(inflowangle, fnormi, Urotor, Utan, rho, twist, polar_alpha, polar_cl, polar_cd)
            
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    vmag2 = Urotor**2 + Utan**2
    lift = 0.5*vmag2*rho*cl*chord
    drag = 0.5*vmag2*rho*cd*chord
    fnorm = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
    ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
    gamma = 0.5*np.sqrt(vmag2)*cl*chord
    return [a, aline, r_R, fnorm, ftan, alpha, inflowangle, twist, chord]


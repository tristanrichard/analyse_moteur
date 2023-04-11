import math as math
import numpy as np
import matplotlib.pyplot as plt

rpm = 3000
s = 1

theta = np.linspace(-180,180,100000)
thetaC = 40
deltaThetaC = 70 
#=== Toutes les données ===#
tau = 10.5     # Vmax/Vmin
D = 0.1        # diametre piston
C = 0.1        # Course
L = 0.15       # Longueur bielle
mpiston = 0.750  # valeur masse piston [kg]
mbielle = 0.600  # valeur masse bielle [kg]  
Gamma = 1.3 
Q=2800000


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """
    Parameters
    ----------
    rpm : type=int
          vitesse de rotation du piston
    s : type=int
        taux de suralimentation
    theta : type=liste numpy
            liste des angles au cours de la rotation
    thetaC : type=int
             angle d'allumage
    deltaThetaC : type=int
                  durée de combustion
    
    Returns 
    -------
    V_output : type=liste numpy
               évolution du volume du cylindre en fonction de l'angle theta
    Q_output : type=liste numpy
               évolution de l'apport de chaleur en fonction de l'angle theta
    F_pied_output : type=liste numpy
                    évolution de la force s'appliquant au pied de la bielle en fonction de l'angle theta
    F_tete_output : type=liste numpy
                    évolution de la force s'appliquant à la tete de la bielle en fonction de l'angle theta
    p_output : type=liste numpy
               évolution de la pression en fonction de l'angle theta
    t : type=int
        epaisseur critique
    
    """
    R = C / 2 # Longueur manivelle
    Beta = L / R
    P_ouput_init=s*100000
    Vc = (2 * R * math.pi * D**2)/4
    thetaRad=np.array(list(map(math.radians,theta)))
    thetaCRad=math.radians(-thetaC)
    deltaThetaCRad=math.radians(deltaThetaC)
    Q_tot = s * (10 ** 5) * Vc / (8.314 * (30 + 273.15))*Q*0.02896 # valeur chaleur emise par fuel par kg de melance admis [J/kg_inlet gas]
    
    
    
    
    def Apport_Q(Q, thetaRad, thetaCRad, deltaThetaCRad):
        """
        Parameters
        ----------
        Q = type=int
            valeur chaleur emise par fuel par kg de melance admis [J/kg_inlet gas]
        thetaRad : type=liste numpy
                liste des angles au cours de la rotation en Rad
        thetaCRad : type=int
                 angle d'allumage en radians
        deltaThetaCRad : type=int
                         durée de combustion en radians
    
        Returns 
        -------
        Ap_Q : type=liste numpy
               évolution de l'apport de chaleur en fonction de l'angle theta
        """
        Ap_Q = np.zeros(len(thetaRad))
        for i in range(len(thetaRad)):
            if thetaCRad <= thetaRad[i] and thetaRad[i] <= thetaCRad + deltaThetaCRad :
                Ap_Q[i]=(Q * (1/2)) *(1- math.cos(math.pi*(thetaRad[i] - thetaCRad)/deltaThetaCRad))
                i=i+1
            else :
                i=i+1 #ap_Q est deja une liste de zéros
        return Ap_Q
    
    
    def Pression(Q, Vc, Beta, gamma, tau, thetaRad, thetaCRad, deltaThetaCRad):
        """
        Parameters
        ----------
        Q = type=int
            valeur chaleur emise par fuel par kg de melance admis [J/kg_inlet gas]
        Vc = type=int
             volume balayé par le piston
        beta = type=int
               rapport entre la longeur de la bielle et la longeur de la manivelle
        gamma = type=int
                constante calorifique dépendante du fluel
        tau = type=int
              rapport entre le volume maximal et minimal du cylindre
        thetaRad : type=liste numpy
                liste des angles au cours de la rotation en Rad
        thetaCRad : type=int
                 angle d'allumage en radians
        deltaThetaCRad : type=int
                         durée de combustion en radians
    
        Returns 
        -------
        Ap_Q : type=liste numpy
               évolution de l'apport de chaleur en fonction de l'angle theta
        """
        V = np.zeros(len(thetaRad))
        for i in range(len(thetaRad)):
            V[i]=((Vc/2)*(1-math.cos(thetaRad[i]) + Beta - math.sqrt(Beta**2 - math.sin(thetaRad[i])**2)) + Vc/(tau-1))     
        
        
        dV = np.zeros(len(thetaRad))
        for i in range(len(thetaRad)):
                dV[i]=(Vc/2 * (math.sin(thetaRad[i]) + (math.sin(thetaRad[i]) * math.cos(thetaRad[i]))/(math.sqrt((Beta*2) - math.sin(thetaRad[i])*2))))
      
        dQ = np.zeros(len(thetaRad))
        for i in range(len(thetaRad)):
            if thetaRad[i] >= thetaCRad and thetaRad[i] <= thetaCRad + deltaThetaCRad :
                dQ[i]=(Q * (1/2) * math.sin( math.pi * (thetaRad[i] - thetaCRad)/deltaThetaCRad) * (math.pi/deltaThetaCRad))
        
        p_output = np.zeros(len(thetaRad))
        for i in range(len(thetaRad)):
            if i>0:
                p_output[i] = p_output[i-1] + (thetaRad[1]-thetaRad[0]) * ((-gamma) * p_output[i - 1] * (dV[i] / V[i]) + ((gamma - 1) * dQ[i] / V[i]))
            else:
                p_output[i] = P_ouput_init + (thetaRad[1]-thetaRad[0]) * ((-gamma) * p_output[i - 1] * (dV[i] / V[i]) + ((gamma - 1) * dQ[i] / V[i]))
        return (p_output, V, dV, dQ)
    
    
    def Force(mpiston, mbielle, thetaRad,P):
        """
        Parameters
        ----------
        mpiston= type=int
                 masse du piston
        mbielle : type=int
                 masse de la bielle
        thetaRad : type=liste numpy
                   liste des angles au cours de la rotation en Rad
        P : type=liste numpy
            pression dans le cylindre en fonction de theta
    
        Returns 
        -------
        F_pied_output : type=liste numpy
                        évolution de la force s'appliquant au pied de la bielle en fonction de l'angle theta
        F_tete_output : type=liste numpy
                        évolution de la force s'appliquant à la tete de la bielle en fonction de l'angle theta
        """
        w=(2*math.pi*rpm)/60
        F_pied_output = np.ones(len(thetaRad))
        for i in range(len(thetaRad)):
            F_pied_output[i]=((math.pi * (D**2) * P[i])/4) - (mpiston * R * (w**2) * math.cos(thetaRad[i]))
        
        F_tete_output=np.ones(len(thetaRad))
        for i in range(len(thetaRad)):
            F_tete_output[i]=(-1) * ((math.pi * (D**2) * P[i])/4) + ((mpiston + mbielle) * R * (w**2) * (math.cos(thetaRad[i])))
            
        return(F_pied_output,F_tete_output)
    
    
    
    def myfunc_t( F_pied_output, F_tete_output):
        """
        Parameters
        ----------
        F_pied_output : type=liste numpy
                        évolution de la force s'appliquant au pied de la bielle en fonction de l'angle theta
        F_tete_output : type=liste numpy
                        évolution de la force s'appliquant à la tete de la bielle en fonction de l'angle theta
    
        Returns 
        -------
        t : type=int
            epaisseur critique
        """
        Fcrit = np.max(np.minimum(F_pied_output, -F_tete_output)) 
        print(Fcrit)
        sigma= 450 * (10**6)
        E = 200 * (10**9) 
        #si flambage dans le plan xx
        Kx=1
        px= np.poly1d([1/Fcrit, 0, -1/(11*sigma), 0, -(Kx*L)**2/(math.pi*math.pi*E*419/12)])
        tx=px.r
        #si flambage dans le plan yy
        Ky=0.5
        py= np.poly1d([1/Fcrit, 0, -1/(11*sigma), 0, -(Ky*L)**2/(math.pi*math.pi*E*131/12)])
        ty=py.r
        txx=[]; tyy=[]
        for elem in tx:
            if np.imag(elem)==0 and elem>0:
                txx.append(elem)
        for elem in ty:
            if np.imag(elem)==0 and elem>0:
                tyy.append(elem)
        return max(abs(min(txx)), abs(min(tyy)))    
     
    
    
    V_output = Pression(Q, Vc, Beta, Gamma, tau, thetaRad, thetaCRad, deltaThetaCRad)[1]
    
    Q_output = Apport_Q(Q_tot, thetaRad, thetaCRad, deltaThetaCRad)
    
    p_output = Pression(Q_tot, Vc, Beta, Gamma, tau, thetaRad, thetaCRad, deltaThetaCRad)[0]
    
    F_pied_output = Force( mpiston, mbielle, thetaRad,p_output)[0]
    F_tete_output = Force( mpiston, mbielle, thetaRad,p_output)[1]
    
    t = myfunc_t( F_pied_output, F_tete_output)
    
    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)

def plot():
    V_output, Q_output, F_pied_output, F_tete_output, p_output, t= myfunc(rpm, s, theta, thetaC, deltaThetaC)
    plt.figure()
    plt.title("pression")
    plt.plot(theta,p_output)
    plt.show()
    
    plt.figure()
    plt.title("volume")
    plt.plot(theta,V_output)
    plt.show()
    
    plt.figure()
    plt.title("Q")
    plt.plot(theta,Q_output)
    plt.show()
    
    plt.figure()
    plt.title("F_pied")
    plt.plot(theta,F_pied_output)
    plt.show()
    
    plt.figure()
    plt.title("F_tete")
    plt.plot(theta,F_tete_output)
    plt.show()
    
    print(t)
    
plot()
    
    
    
    
    
    
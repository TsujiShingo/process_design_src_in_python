import numpy as np
from scipy import integrate #熱容量の積分

R = 8.314462 #[J K^-1 mol^-1]　# 気体定数
Compounds = ("co","co2","h2","h2o","meoh","n2") #登場する化合物一覧


class Reactor:
    def __init__(self,Pressure=82,catalyst=1100,D=44.5*(10**-3),L=7260*(10**-3),n=4801): # nは個数
        print("Let's costruct Reactor !")
        self.Pressure = Pressure                # [bar] # 圧力
        self.catalyst = catalyst                # [kg m^-3] # 触媒の充填密度
        self.D        = D                       # [m] # 反応器断面の直径
        self.L        = L                       # [m] #反応器長さ
        self.S        = np.pi * (D**2)/4        # [m^2] # 反応器断面積
        self.vmax     = self.S*L                     # [m^3] # 反応器の大きさ
        print("Pressure : ", self.Pressure," bar\n",
            "catalyst : ", self.catalyst, " kg / m^3\n",
            "D : ", self.D, " m\n",
            "L : ", self.L, " m\n",
            "S : ", self.S, " m^2 \n",
            "vmax : ",self.vmax, "m^3 \n" 
        )

class Flow:
    def __init__(self,f_co=0,f_co2=0,f_h2=0,f_h2o=0,f_meoh=0,f_n2=0):
        self.Molar_flow =dict() 
        self.Molar_flow["co"]   = f_co
        self.Molar_flow["co2"]  = f_co2
        self.Molar_flow["h2"]   = f_h2
        self.Molar_flow["h2o"]  = f_h2o
        self.Molar_flow["meoh"] = f_meoh
        self.Molar_flow["n2"]   = f_n2

    def y(self,material):#モル分率を返す
        sum_Molar_flow = 0
        for compound in Compounds:
            sum_Molar_flow = sum_Molar_flow + self.Molar_flow[compound]
        if sum_Molar_flow == 0:
            print("sum_Molar_flow is 0 !")
            return 0
        return self.Molar_flow[material]/sum_Molar_flow


## Reactions 
def K_MeOH(T): # [mol kg_cat^-1 s^-1 bar^-2]
    return 1.07*np.exp(36696/R/T)

def K_RWGS(T): # [mol kg_cat^-1 s^-1 bar^-1]
    return 1.22*(10**10)*np.exp(-94765/R/T)

def K_a(T): # [bar^-0.5]
    return 0.499*np.exp(17197/R/T)

def K_b(T): # [bar^-1]
    return 6.62*(10**-11)*np.exp(124119/R/T)

def K_c(T):#[-]
    return 3453.38

def K_eq1(T): # [bar^-2]
    return 10**(3066/T-10.592)

def K_eq2(T): #[-]
    return 10**(-2073/T+2.029)

def reaction_rate_1(T,reactor,flow): # [ mol kg_cat^-1 s^-1]
    p_co = reactor.Pressure*flow.y("co")
    p_co2 = reactor.Pressure*flow.y("co2")
    p_h2 = reactor.Pressure*flow.y("h2")
    p_h2o = reactor.Pressure*flow.y("h2o")
    p_meoh = reactor.Pressure*flow.y("meoh")
    return K_MeOH(T)*p_co2*p_h2*(1-p_meoh*p_h2o/(K_eq1(T)*(p_h2**3)*p_co2))/ \
            ((1+K_c(T)*p_h2o/p_h2+K_a(T)*(p_h2**0.5)+K_b(T)*p_h2o)**3)

def reaction_rate_2(T,reactor,flow): # [ mol kg_cat^-1 s^-1]
    p_co = reactor.Pressure*flow.y("co")
    p_co2 = reactor.Pressure*flow.y("co2")
    p_h2 = reactor.Pressure*flow.y("h2")
    p_h2o = reactor.Pressure*flow.y("h2o")
    p_meoh = reactor.Pressure*flow.y("meoh")        
    return K_RWGS(T)*p_co2*(1-p_co*p_h2o/(K_eq2(T)*p_h2*p_co2))/ \
            (1+K_c(T)*p_h2o/p_h2+K_a(T)*(p_h2**0.5)+K_b(T)*p_h2o)

## Heat
dH1_298K = -49.5*(10**3) #反応1での標準反応熱 [J mol^-1]
dH2_298K = 41.2*(10**3)  #反応2での標準反応熱 [J mol^-1]

def Cp_co(T): #[J K^-1 mol^-1]
    return 30.871+(-0.0125)*T+(2.79*(10**-5))*(T**2)+(-1.27*(10**-8))*(T**3)

def Cp_co2(T): #[J K^-1 mol^-1]
    return 19.796+(7.34*(10**-2))*T+(-5.60*(10**-5))*(T**2)+(1.72*(10**-8))*(T**3)

def Cp_h2(T): #[J K^-1 mol^-1]
    return 27.144+(-9.27*(10**-3))*T+(-1.38*(10**-5))*(T**2)+(7.65*(10**-9))*(T**3)

def Cp_h2o(T):#[J K^-1 mol^-1]
    return 32.244+(1.92*(10**-3))*T+(1.06*(10**-5))*(T**2)+(-3.60*(10**-9))*(T**3)

def Cp_meoh(T): #[J K^-1 mol^-1]
    return 21.153+(7.09*(10**-2))*T+(2.59*(10**-5))*(T**2)+(-2.85*(10**-8))*(T**3)

def Cp_n2(T):
    return 31.151+(-1.357*(10**-2))*T+(2.680*(10**-5))*(T**2)+(-1.168*(10**-8))*(t**3)

#温度変化による反応熱(反応エンタルピー)の計算
def Heat_of_reaction(T,Xi1,Xi2):#T K基準でr1の反応がxi1[mol/s],r2の反応がxi2[mol/s]だけ進んだときの反応熱 [J/s]
    H_co,err_co=integrate.quad(Cp_co,298,T)
    H_co2,err_co2=integrate.quad(Cp_co2,298,T)
    H_h2,err_h2=integrate.quad(Cp_h2,298,T)
    H_h2o,err_h2o=integrate.quad(Cp_h2o,298,T)
    H_meoh,err_meoh=integrate.quad(Cp_meoh,298,T)
    return (dH1_298K+H_meoh+H_h2o-H_co2-3*H_h2)*Xi1+(dH2_298K+H_co+H_h2o-H_co2-H_h2)*Xi2

#298K基準での各組成・温度でのエンタルピー計算
def Sum_FH(T,flow):
    H_co,err_co=integrate.quad(Cp_co,298,T)
    H_co2,err_co2=integrate.quad(Cp_co2,298,T)
    H_h2,err_h2=integrate.quad(Cp_h2,298,T)
    H_h2o,err_h2o=integrate.quad(Cp_h2o,298,T)
    H_meoh,err_meoh=integrate.quad(Cp_meoh,298,T)
    H_n2,err_meoh=integrate.quad(Cp_n2,298,T)
    return flow.Molar_flow["co"]*H_co+flow.Molar_flow["co2"]*H_co2+\
           flow.Molar_flow["h2"]*H_h2+flow.Molar_flow["h2o"]*H_h2o+\
           flow.Molar_flow["meoh"]*H_meoh+flow.Molar_flow["n2"]*H_n2


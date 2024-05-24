import numpy as np
from scipy import integrate #熱容量の積分

R = 8.314462 #[J K^-1 mol^-1]　# 気体定数
Compounds = ("co","co2","h2","h2o","meoh","n2") #登場する化合物一覧 必ず小文字で
stoichiometry_matrix = { #量論係数 (反応1,反応2)
    "co"    :( 0, 1),
    "co2"   :(-1,-1),
    "h2"    :(-3,-1),
    "h2o"   :( 1, 1),
    "meoh"  :( 1, 0),
    "n2"    :( 0, 0)
}


class Reactor:
    def __init__(self,Pressure=82,catalyst=1100,D=44.5*(10**-3),L=7260*(10**-3)): # nは個数
        print("Let's costruct Reactor !")
        self.Pressure = Pressure                # [bar] # 圧力
        self.catalyst = catalyst                # [kg m^-3] # 触媒の充填密度
        self.D        = D                       # [m] # 反応器断面の直径
        self.L        = L                       # [m] #反応器長さ
        self.S        = np.pi * (D**2)/4        # [m^2] # 反応器断面積
        self.vmax     = L*self.S                     # [m^3] # 反応器の大きさ
        print("Pressure : ", self.Pressure," bar\n",
            "catalyst : ", self.catalyst, " kg / m^3\n",
            "D : ", self.D, " m\n",
            "L : ", self.L, " m\n",
            "S : ", self.S, " m^2 \n",
            "vmax : ",self.vmax, "m^3 \n" 
        )

class Flow: # [mol/s] #辞書型で初期化可能
    def __init__(self, molar_flows = None):
        default_flows = {compound : 0 for compound in Compounds} # 0 初期化
        if molar_flows is None:
            molar_flows = default_flows
        else: # 登録外の化合物を見つけたらエラー
            for key in molar_flows:
                if key not in default_flows:
                    raise ValueError(f"Unknown key in molar_flows: {key}")
            # 2つ辞書をマージ,同じものがあれば右側にある初期値に上書き
            self.Molar_flow = {**default_flows, **molar_flows}
        # self.Molar_flowを設定
        self.Molar_flow = {compound: molar_flows.get(compound, 0) for compound in Compounds}

    def total_flow(self): # 全モル流量
        return sum(self.Molar_flow.values())

    def y(self,compound): # モル分率
        if compound not in self.Molar_flow:
            raise ValueError(f"Compound {compound} not found in flow.")
        if self.total_flow() <= 0:
            raise ValueError(f"total_flow not Positive.")
        else :
            return self.Molar_flow[compound] / self.total_flow()

class ReactionRates:
    @staticmethod
    def K_MeOH(T): # [mol kg_cat^-1 s^-1 bar^-2]
        return 1.07*np.exp(36696/R/T)

    @staticmethod
    def K_RWGS(T): # [mol kg_cat^-1 s^-1 bar^-1]
        return 1.22*(10**10)*np.exp(-94765/R/T)

    @staticmethod
    def K_a(T): # [bar^-0.5]
        return 0.499*np.exp(17197/R/T)

    @staticmethod
    def K_b(T): # [bar^-1]
        return 6.62*(10**-11)*np.exp(124119/R/T)

    @staticmethod
    def K_c(T):#[-]
        return 3453.38

    @staticmethod
    def K_eq1(T): # [bar^-2]
        return 10**(3066/T-10.592)

    @staticmethod
    def K_eq2(T): #[-]
        return 10**(-2073/T+2.029)

    @staticmethod
    def reaction_rate_1(T,reactor,flow): # [ mol kg_cat^-1 s^-1]
        p = dict()
        for compound in Compounds:
            p[compound] = reactor.Pressure*flow.y(compound)    
        return ReactionRates.K_MeOH(T)*p["co2"]*p["h2"]*(1-p["meoh"]*p["h2o"]/(ReactionRates.K_eq1(T)*(p["h2"]**3)*p["co2"]))/ \
        ((1+ReactionRates.K_c(T)*p["h2o"]/p["h2"]+ReactionRates.K_a(T)*(p["h2"]**0.5)+ReactionRates.K_b(T)*p["h2o"])**3)
    
    @staticmethod
    def reaction_rate_2(T,reactor,flow): # [ mol kg_cat^-1 s^-1]
        p = dict()
        for compound in Compounds:
            p[compound] = reactor.Pressure*flow.y(compound)    
        return ReactionRates.K_RWGS(T)*p["co2"]*(1-p["co"]*p["h2o"]/(ReactionRates.K_eq2(T)*p["h2"]*p["co2"]))/ \
                (1+ReactionRates.K_c(T)*p["h2o"]/p["h2"]+ReactionRates.K_a(T)*(p["h2"]**0.5)+ReactionRates.K_b(T)*p["h2o"])

class CpFunctions_T3: #[J K^-1 mol^-1] 
    # a + b*T + c*T^2 + d*T^3型で書けることを利用する
    Cp_coefficient = { #(a,b,c,d)形式で保存
        "co"    :(30.871,        -0.0125, 2.79*(10**-5), -1.27*(10**-8)),
        "co2"   :(19.796,  7.34*(10**-2),-5.60*(10**-5),  1.72*(10**-8)),
        "h2"    :(27.144, -9.27*(10**-3),-1.38*(10**-5),  7.65*(10**-9)),
        "h2o"   :(32.244,  1.92*(10**-3), 1.06*(10**-5), -3.60*(10**-9)),
        "meoh"  :(21.153,  7.09*(10**-2), 2.59*(10**-5), -2.85*(10**-8)),
        "n2"    :(31.151,-1.357*(10**-2),2.680*(10**-5),-1.168*(10**-8))
    }

    @staticmethod
    def CpFunc_T3(compound,T):
        coeffs = CpFunctions_T3.Cp_coefficient[compound]
        return np.inner(coeffs,(1,T,T**2,T**3))

class EnthalpyCalculations:#温度変化による反応熱(反応エンタルピー)の計算
    @staticmethod
    def Heat_of_reaction(T,Xi1,Xi2):#T K基準でr1の反応がXi1[mol/s],r2の反応がXi2[mol/s]だけ進んだときの反応熱 [J/s]
        dH1_298K = -49.5*(10**3) #反応1での標準反応熱 [J mol^-1]
        dH2_298K = 41.2*(10**3)  #反応2での標準反応熱 [J mol^-1]
        Cp = CpFunctions_T3()
        H=dict()
        for compound in Compounds:
            H[compound],_ = integrate.quad(lambda T:Cp.CpFunc_T3(compound,T),298,T) 
        delta_H1 = sum(stoichiometry_matrix[compound][0] * H[compound] for compound in Compounds)
        delta_H2 = sum(stoichiometry_matrix[compound][1] * H[compound] for compound in Compounds)
        return (dH1_298K + delta_H1)*Xi1 + (dH2_298K + delta_H2)*Xi2

    @staticmethod
    #298K基準での各組成・温度でのエンタルピー計算
    def Sum_FH(T,flow):
        Cp = CpFunctions_T3()
        H=dict()
        for compound in Compounds:
            H[compound],_ = integrate.quad(lambda T:Cp.CpFunc_T3(compound,T),298,T)
        return sum(flow.Molar_flow[compound] * H[compound] for compound in Compounds)

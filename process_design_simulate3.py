import numpy as np
from scipy.optimize import fsolve #熱バランス解くときに使った
import process_design_module3 as pdm #自作
from process_design_module3 import ReactionRates
from process_design_module3 import EnthalpyCalculations

# 反応器の作成
# def create_reactor():
#     a0 = 3 # [m s^-1] # 空塔速度
#     u0 = a0*reactor.S #[m^3/s]　#体積速度
def progress_debug(n,n_step):
    print(f"\r {n}/{n_step}",end="")

def simulate_no_hysys_isothermal(flow_in,T_in,reactor,dv,n_step,progress_callback = progress_debug):
    flows = list() 
    flows.append(flow_in) #初期値設定
    Temperature = list()
    Temperature.append(T_in)
    Volume = list()
    Volume.append(0)
    reaction = ReactionRates()
    for n in range(n_step-1): #シミュレーション開始
        flows.append(pdm.Flow())
        Temperature.append(T_in)
        Volume.append(round(dv*(n+1),2))
        xi1 = reaction.reaction_rate_1(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
        xi2 = reaction.reaction_rate_2(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
        #モル流量の更新 F_out = F_in + (a*r1+b*r2)dv
        flows[n+1].Molar_flow["co"]   = flows[n].Molar_flow["co"]   + ( 0*xi1) + ( 1*xi2)
        flows[n+1].Molar_flow["co2"]  = flows[n].Molar_flow["co2"]  + (-1*xi1) + (-1*xi2)
        flows[n+1].Molar_flow["h2"]   = flows[n].Molar_flow["h2"]   + (-3*xi1) + (-1*xi2)
        flows[n+1].Molar_flow["h2o"]  = flows[n].Molar_flow["h2o"]  + ( 1*xi1) + ( 1*xi2)
        flows[n+1].Molar_flow["meoh"] = flows[n].Molar_flow["meoh"] + ( 1*xi1) + ( 0*xi2)
        flows[n+1].Molar_flow["n2"]   = flows[n].Molar_flow["n2"] #不活性
        #熱の変更なし
        Temperature[n+1]=Temperature[n]
        progress_callback(n+1,n_step)
    print("simulate_no_hysys_isothermal finish")
    return (Volume,Temperature,flows)

def simulate_no_hysys_adiabatic(flow_in,T_in,reactor,dv,n_step,progress_callback = progress_debug):
    flows = list() 
    flows.append(flow_in) #初期値設定
    Temperature = list()
    Temperature.append(T_in)
    Volume = list()
    Volume.append(0)
    reaction = ReactionRates()
    EnthalpyCal = EnthalpyCalculations()
    for n in range(n_step-1): #シミュレーション開始
        flows.append(pdm.Flow())
        Temperature.append(T_in)
        Volume.append(round(dv*(n+1),2))
        xi1 = reaction.reaction_rate_1(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
        xi2 = reaction.reaction_rate_2(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
        #モル流量の更新 F_out = F_in + (a*r1+b*r2)dv
        flows[n+1].Molar_flow["co"]   = flows[n].Molar_flow["co"]   + ( 0*xi1) + ( 1*xi2)
        flows[n+1].Molar_flow["co2"]  = flows[n].Molar_flow["co2"]  + (-1*xi1) + (-1*xi2)
        flows[n+1].Molar_flow["h2"]   = flows[n].Molar_flow["h2"]   + (-3*xi1) + (-1*xi2)
        flows[n+1].Molar_flow["h2o"]  = flows[n].Molar_flow["h2o"]  + ( 1*xi1) + ( 1*xi2)
        flows[n+1].Molar_flow["meoh"] = flows[n].Molar_flow["meoh"] + ( 1*xi1) + ( 0*xi2)
        flows[n+1].Molar_flow["n2"]   = flows[n].Molar_flow["n2"] #不活性
        #熱の変更
        def heat_balance(T):
            return EnthalpyCal.Sum_FH(T,flows[n+1])- \
                    EnthalpyCal.Sum_FH(Temperature[n],flows[n+1])+ \
                    EnthalpyCal.Heat_of_reaction(Temperature[n],xi1,xi2)
        Temperature[n+1],*_ =fsolve(heat_balance,Temperature[n]) #heat_balanece(T)=0 となるTを求める
        progress_callback(n+1,n_step)
    print("simulate_no_hysys_adiabatic finish")
    return (Volume,Temperature,flows)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import pandas
    reactor = pdm.Reactor(Pressure=30,D=2.333599687,L=7.26) #標準では 論文と同じ条件で作成
    dv=0.01 
    n_step = int(reactor.vmax/dv)
    flow_in = pdm.Flow({"h2":29520.0/3.6,"n2":997.9/3.6,"co":3017.8/3.6,"co2":10800.0/3.6,"meoh":674.4/3.6,"h2o":172.9/3.6}) #辞書形式で初期値与えて
    T_in = 485 
    Volume,Temperature,flows=simulate_no_hysys_adiabatic(flow_in,T_in,reactor,dv,n_step)
    data = list() #微小体積毎の温度・モル流量をcsvfileで保存
    for n in range(n_step): #各体積での処理
        data.append({"V[m^3]":Volume[n],"T[K]":Temperature[n],"f_co[mol/s]":flows[n].Molar_flow["co"],"f_co2[mol/s]":flows[n].Molar_flow["co2"],"f_h2[mol/s]":flows[n].Molar_flow["h2"],"f_h2o[mol/s]":flows[n].Molar_flow["h2o"],"f_meoh[mol/s]":flows[n].Molar_flow["meoh"],"f_n2[mol/s]":flows[n].Molar_flow["n2"]})
    data_csv = pandas.DataFrame(data)
    data_csv.plot(grid=True,x="V[m^3]",y="T[K]")
    plt.savefig("温度変化.png")
    plt.show()

    data_csv[["V[m^3]","f_co[mol/s]","f_co2[mol/s]","f_h2[mol/s]","f_h2o[mol/s]","f_meoh[mol/s]","f_n2[mol/s]"]].plot(
    grid=True,x="V[m^3]")
    plt.savefig("流量変化(全体).png")
    plt.show()

    data_csv[["V[m^3]","f_co[mol/s]","f_co2[mol/s]","f_h2[mol/s]","f_h2o[mol/s]","f_meoh[mol/s]","f_n2[mol/s]"]].plot(
        subplots = True,grid=True,x="V[m^3]")
    plt.savefig("流量変化(個別).png")
    plt.show()

    conditions = pandas.DataFrame([{"Pressure[bar]":reactor.Pressure,"catalyst[kg/m^3]" : reactor.catalyst,\
        "D[m]":reactor.D,"L[m]":reactor.L,"S[m^2]":reactor.S,"vmax [m^3]":reactor.vmax,\
        "dv[m^3]":dv,"n_step":n_step}]).T

    conditonsname = 'conditions.txt'
    conditions.to_csv(conditonsname)
    datafile = "datafile.csv"
    data_csv.to_csv(datafile,index=False)



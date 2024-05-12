#%matplotlib inline #jupyter notebookでは必要
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate #熱容量の積分
from scipy.optimize import fsolve #熱バランス解くときに使った
import sys
import pandas #csvファイルの入出力系
import datetime
import os
import process_design_module as pdm #自作

# 反応器の作成
reactor = pdm.Reactor() #標準では 論文と同じ条件で作成
reactor.S = reactor.S * 4801.0  #4801個の反応器の並列分だけ断面積を変化
reactor.vmax = reactor.vmax * 4801.0 
a0 = 3 # [m s^-1] # 空塔速度
u0 = a0*reactor.S #[m^3/s]　#体積速度

# シミュレーションのステップ数の決定
dv = 0.01 #[m^3] # 微小体積
n_step = int(reactor.vmax/dv) #ステップ数
# 微小体積の流量の情報を保存
flows=list() 
for n in range(n_step):
    flows.append(pdm.Flow())

Temperature=np.full(n_step,500.0)

# 0.実行結果を保存するディレクトリの作成
now = datetime.datetime.now()
foldername= "output/" + now.strftime('%Y%m%d_%H%M%S')+"/"
os.makedirs(foldername, exist_ok=True) #保存先の作成

# 1.入口流量の初期化
Ft0 =  reactor.Pressure*(10**5)*u0/(pdm.R*Temperature[0]) # PV=nRT -> F=Pu/RT
y0_co, y0_co2, y0_h2, y0_h2o, y0_meoh = (0.017, 0.173, 0.803, 0.001, 0.005)

flows[0].Moler_flow["co"]   = Ft0 * y0_co
flows[0].Moler_flow["co2"]  = Ft0 * y0_co2
flows[0].Moler_flow["h2"]   = Ft0 * y0_h2
flows[0].Moler_flow["h2o"]  = Ft0 * y0_h2o
flows[0].Moler_flow["meoh"] = Ft0 * y0_meoh


# mid-output
print("入口組成：T = {} K,  bar,\n F_co = {} mol/s, F_co2 = {} mol/s, F_h2 ={} mol/s,\n F_h2o = {} mol/s, F_meoh = {} mol/s\n".format(
    Temperature[0],flows[0].Moler_flow["co"],flows[0].Moler_flow["co2"],flows[0].Moler_flow["h2"],flows[0].Moler_flow["h2o"],flows[0].Moler_flow["meoh"]))

# 2. 管内のシミュレーション
for n in range(n_step-1):
    #反応速度式の計算
    xi1 = pdm.reaction_rate_1(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
    xi2 = pdm.reaction_rate_2(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
    #モル流量の更新 F_out = F_in + (a*r1+b*r2)dv
    flows[n+1].Moler_flow["co"]   = flows[n].Moler_flow["co"]   + ( 0*xi1) + ( 1*xi2)
    flows[n+1].Moler_flow["co2"]  = flows[n].Moler_flow["co2"]  + (-1*xi1) + (-1*xi2)
    flows[n+1].Moler_flow["h2"]   = flows[n].Moler_flow["h2"]   + (-3*xi1) + (-1*xi2)
    flows[n+1].Moler_flow["h2o"]  = flows[n].Moler_flow["h2o"]  + ( 1*xi1) + ( 1*xi2)
    flows[n+1].Moler_flow["meoh"] = flows[n].Moler_flow["meoh"] + ( 1*xi1) + ( 0*xi2)
    #熱の変更
    def heat_balance(T):
        return pdm.Sum_FH(T,flows[n+1])- \
                pdm.Sum_FH(Temperature[n],flows[n+1])+ \
                pdm.Heat_of_reaction(Temperature[n],xi1,xi2)
#   temperature[n+1]=temperature[n] #等温の場合
    Temperature[n+1],*_ =fsolve(heat_balance,Temperature[n]) #heat_balanece(T)=0 となるTを求める
    # mid-output
    sys.stdout.write("\r simulation:step: {} / {}, T = {} K".format(n+1,n_step,Temperature[n+1]))
    sys.stdout.flush()

# mid-output 
print("\n出口組成：T = {} K,\n F_co = {} mol/s, F_co2 = {} mol/s, F_h2 ={} mol/s,\n F_h2o = {} mol/s, F_meoh = {} mol/s\n".format(
    Temperature[-1],flows[-1].Moler_flow["co"],flows[-1].Moler_flow["co2"],flows[-1].Moler_flow["h2"],flows[-1].Moler_flow["h2o"],flows[-1].Moler_flow["meoh"]))

# #3.シミュレーション結果の出力
conditions = pandas.DataFrame([{"Pressure[bar]":reactor.Pressure,"catalyst[kg/m^3]" : reactor.catalyst,\
    "D[m]":reactor.D,"L[m]":reactor.L,"S[m^2]":reactor.S,"vmax [m^3]":reactor.vmax,\
    "dv[m^3]":dv,"Superficial velocity (inlet) a0[m/s]":a0,"n_step":n_step}]).T

conditonsname = foldername+'conditions.txt'
conditions.to_csv(conditonsname)
print(conditions)

data = list() #微小体積毎の温度・モル流量をcsvfileで保存
for n in range(n_step): #各体積での処理
    data.append({"V[m^3]":dv*n,"T[K]":Temperature[n],"f_co[mol/s]":flows[n].Moler_flow["co"],"f_co2[mol/s]":flows[n].Moler_flow["co2"],"f_h2[mol/s]":flows[n].Moler_flow["h2"],"f_h2o[mol/s]":flows[n].Moler_flow["h2o"],"f_meoh[mol/s]":flows[n].Moler_flow["meoh"]})

data_csv = pandas.DataFrame(data)
datafile = foldername+ "datafile.csv"
data_csv.to_csv(datafile,index=False)
print(data_csv)

#モル流量のグラフ
data_csv[["V[m^3]","f_co[mol/s]","f_co2[mol/s]","f_h2[mol/s]","f_h2o[mol/s]","f_meoh[mol/s]"]].plot(
    grid=True,x="V[m^3]")
plt.savefig(foldername+"Moler_flow.png")
plt.show()

data_csv[["V[m^3]","f_co[mol/s]","f_co2[mol/s]","f_h2[mol/s]","f_h2o[mol/s]","f_meoh[mol/s]"]].plot(
    subplots = True,grid=True,x="V[m^3]")
plt.savefig(foldername+"Moler_flow(detail).png")
plt.show()

#温度分布
data_csv.plot(grid=True,x="V[m^3]",y="T[K]")
plt.savefig(foldername+"Temperature.png")
plt.show()

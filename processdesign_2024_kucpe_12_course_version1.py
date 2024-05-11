#%matplotlib inline #jupyter notebookでは必要
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate #熱容量の積分
from scipy.optimize import fsolve #熱バランス解くときに使った
import sys
import pandas #csvファイルの入出力系
import datetime
import os

R = 8.314462 #[J K^-1 mol^-1]　# 気体定数
catalyst = 1100 #[kg m^-3]　# 触媒の充填密度

D = 2.8 #[m] # 反応器断面の直径
S = np.pi * (D**2) /4 #[m^2]反応器断面積
L = 2*D #[m] #反応器長さ
vmax = S*L # [m^3] # 反応器の大きさ
a = 3 # [m s^-1] #
dv = 0.01 #[m^3] # 微小体積
u = a*S #[m^3/s]　#空塔速度

n_step = int(vmax/dv) #ステップ数

temperature = np.ones(n_step)# [K] #温度
pressure    = np.zeros(n_step) # [bar](1 bar=10^5 Par) #圧力 

#モル流量
f_co  = np.zeros(n_step) # [mol m^-3 s^-1]
f_co2 = np.zeros(n_step) # [mol m^-3 s^-1]
f_h2  = np.zeros(n_step) # [mol m^-3 s^-1]
f_h2o = np.zeros(n_step) # [mol m^-3 s^-1]
f_meoh= np.zeros(n_step) # [mol m^-3 s^-1]

#分圧 [bar]
p_co  = np.zeros(n_step) # [bar]
p_co2 = np.zeros(n_step) # [bar]
p_h2  = np.zeros(n_step) # [bar]
p_h2o = np.zeros(n_step) # [bar]
p_meoh= np.zeros(n_step) # [bar]

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

def reaction_rate_1(T,V): # [ mol kg_cat^-1 s^-1]
    return K_MeOH(T)*p_co2[V]*p_h2[V]*(1-p_meoh[V]*p_h2o[V]/(K_eq1(T)*(p_h2[V]**3)*p_co2[V]))/ \
            ((1+K_c(T)*p_h2o[V]/p_h2[V]+K_a(T)*(p_h2[V]**0.5)+K_b(T)*p_h2o[V])**3)

def reaction_rate_2(T,V): # [ mol kg_cat^-1 s^-1]
    return K_RWGS(T)*p_co2[V]*(1-p_co[V]*p_h2o[V]/(K_eq2(T)*p_h2[V]*p_co2[V]))/ \
            (1+K_c(T)*p_h2o[V]/p_h2[V]+K_a(T)*(p_h2[V]**0.5)+K_b(T)*p_h2o[V])

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

#温度変化による反応熱(反応エンタルピー)の計算
def Heat_of_reaction(T,Xi1,Xi2):#T K基準でr1の反応がxi1[mol/s],r2の反応がxi2[mol/s]だけ進んだときの反応熱 [J/s]
    H_co,err_co=integrate.quad(Cp_co,298,T)
    H_co2,err_co2=integrate.quad(Cp_co2,298,T)
    H_h2,err_h2=integrate.quad(Cp_h2,298,T)
    H_h2o,err_h2o=integrate.quad(Cp_h2o,298,T)
    H_meoh,err_meoh=integrate.quad(Cp_meoh,298,T)
    return (dH1_298K+H_meoh+H_h2o-H_co2-3*H_h2)*Xi1+(dH2_298K+H_co+H_h2o-H_co2-H_h2)*Xi2

#298K基準での各組成・温度でのエンタルピー計算
def Sum_FH(T,F_co,F_co2,F_h2,F_h2o,F_meoh):
    H_co,err_co=integrate.quad(Cp_co,298,T)
    H_co2,err_co2=integrate.quad(Cp_co2,298,T)
    H_h2,err_h2=integrate.quad(Cp_h2,298,T)
    H_h2o,err_h2o=integrate.quad(Cp_h2o,298,T)
    H_meoh,err_meoh=integrate.quad(Cp_meoh,298,T)
    return F_co*H_co+F_co2*H_co2+F_h2*H_h2+F_h2o*H_h2o+F_meoh*H_meoh

def initialize_molar_flow_rate_and_pressure_byPT(P0,T0,y_co,y_co2,y_h2,y_h2o,y_meoh):
    pressure[0] = P0 #全圧
    temperature[0]=T0
    p_co[0]   = P0*y_co
    p_co2[0]  = P0*y_co2
    p_h2[0]   = P0*y_h2
    p_h2o[0]  = P0*y_h2o
    p_meoh[0] = P0*y_meoh
    f_total = P0*(10**5)*u/(R*T0) # PV=nRT の両辺を単位時間で割って-> Pu=FRT
    f_co[0]   = f_total*y_co # [mol m^-3 s^-1]
    f_co2[0]  = f_total*y_co2 # [mol m^-3 s^-1]
    f_h2[0]   = f_total*y_h2 #[mol m^-3 s^-1]
    f_h2o[0]  = f_total*y_h2o #[mol m^-3 s^-1]
    f_meoh[0] = f_total*y_meoh #[mol m^-3 s^-1]

def setting_pressure(n):
    pressure[n]=(f_co[n]+f_co2[n]+f_h2[n]+f_h2o[n]+f_meoh[n])/u*R*temperature[n]/(10**5) #全圧
    p_co[n]   = f_co[n]/u*R*temperature[n]/(10**5)
    p_co2[n]  = f_co2[n]/u*R*temperature[n]/(10**5)
    p_h2[n]   = f_h2[n]/u*R*temperature[n]/(10**5)
    p_h2o[n]  = f_h2o[n]/u*R*temperature[n]/(10**5)
    p_meoh[n] = f_meoh[n]/u*R*temperature[n]/(10**5)


# 0.実行結果を保存するディレクトリの作成
now = datetime.datetime.now()
foldername= "output/" + now.strftime('%Y%m%d_%H%M%S')+"/"
# os.makedirs(foldername, exist_ok=True) #保存先の作成

# 1.管型反応器の入口条件を設定,出力
initialize_molar_flow_rate_and_pressure_byPT(50,550,0.01,0.48,0.48,0.01,0.02)
# mid-output
print("入口組成：T = {} K, P = {} bar,\n F_co = {} mol/s, F_co2 = {} mol/s, F_h2 ={} mol/s,\n F_h2o = {} mol/s, F_meoh = {} mol/s\n".format(
    temperature[0],pressure[0],f_co[0],f_co2[0],f_h2[0],f_h2o[0],f_meoh[0]))

# 2. 管内のシミュレーション
for n in range(n_step-1):
    #反応速度式の計算
    r1 = reaction_rate_1(temperature[n],n)
    r2 = reaction_rate_2(temperature[n],n)
    xi1 = r1*catalyst*dv
    xi2 = r2*catalyst*dv
    #モル流量の更新 F_out = F_in + (a*r1+b*r2)dv
    f_co[n+1]   = f_co[n]   + ( 0*xi1) + ( 1*xi2)
    f_co2[n+1]  = f_co2[n]  + (-1*xi1) + (-1*xi2)
    f_h2[n+1]   = f_h2[n]   + (-3*xi1) + (-1*xi2)
    f_h2o[n+1]  = f_h2o[n]  + ( 1*xi1) + ( 1*xi2)
    f_meoh[n+1] = f_meoh[n] + ( 1*xi1) + ( 0*xi2)
    #熱の変更
    def heat_balance(T):
        return Sum_FH(T,f_co[n+1],f_co2[n+1],f_h2[n+1],f_h2o[n+1],f_meoh[n+1])- \
                Sum_FH(temperature[n],f_co[n+1],f_co2[n+1],f_h2[n+1],f_h2o[n+1],f_meoh[n+1])+ \
                Heat_of_reaction(temperature[n],xi1,xi2)
#   temperature[n+1]=temperature[n] #等温の場合
    temperature[n+1],*_ =fsolve(heat_balance,temperature[n]) #heat_balanece(T)=0 となるTを求める

    #圧力の更新 F_out 側
    setting_pressure(n+1)
    # mid-output
    sys.stdout.write("\r simulation:step: {} / {}, T = {} K,P={} bar".format(n+1,n_step,temperature[n+1],pressure[n+1]))
    sys.stdout.flush()

# mid-output
print("\n\n出口組成：T = {} K, P = {} bar,\n F_co = {} mol/s, F_co2 = {} mol/s, F_h2 ={} mol/s,\n F_h2o = {} mol/s, F_meoh = {} mol/s\n".format(
    temperature[-1],pressure[-1],f_co[-1],f_co2[-1],f_h2[-1],f_h2o[-1],f_meoh[-1]))

#3.シミュレーション結果の出力
conditions = pandas.DataFrame([{"catalyst [kg m^-3]" : catalyst,"vmax [m^3]":vmax,"dv[m^3]":dv,"n_step":n_step,"Superficial velocity u[m^3/s]":u}])
conditonsname = foldername+'conditions.txt'
# conditions.T.to_csv(conditonsname, index=True,)
print(conditions)

data = list() #各モル流量とかをcsvfileで残す
for n in range(n_step): #各体積での処理
    data.append({"step":n,"V[m^3]":dv*n,"T[K]":temperature[n],"P[bar]":pressure[n],"f_co[mol/s]":f_co[n],"f_co2[mol/s]":f_co2[n],"f_h2[mol/s]":f_h2[n],"f_h2o[mol/s]":f_h2o[n],"f_meoh[mol/s]":f_meoh[n]})

data_csv = pandas.DataFrame(data)
datafile = foldername+ "断熱非等温反応器_試作_ver3_実行結果" + '.csv'
# data_csv.to_csv(datafile)

print(data_csv)

#グラフ
#モル流量
fig,ax = plt.subplots()
ax.set_xlabel('Volume in reactor [m^3]')
ax.set_ylabel('Mass flow rate [mol/s]')
plt.grid()
plt_v = np.linspace(0,vmax,n_step)
plt.plot(plt_v,f_co,label="CO")
plt.plot(plt_v,f_co2,label="CO2")
plt.plot(plt_v,f_h2,label="H2")
plt.plot(plt_v,f_h2o,label="H2O")
plt.plot(plt_v,f_meoh,label="MeOH")
plt.legend()
# plt.savefig(foldername+"Mass_flow_rate.png")   
plt.show()

#圧力
fig2,ax2 = plt.subplots()
ax2.set_xlabel('Volume in reactor [m^3]')
ax2.set_ylabel('Pressure [bar]')
plt.grid()
plt.plot(plt_v,pressure,label="total")
plt.plot(plt_v,p_co,label="CO")
plt.plot(plt_v,p_co2,label="CO2")
plt.plot(plt_v,p_h2,label="H2")
plt.plot(plt_v,p_h2o,label="H2O")
plt.plot(plt_v,p_meoh,label="MeOH")
plt.legend()
# plt.savefig(foldername+"Pressure.png")   
plt.show()

#温度
fig3,ax3 = plt.subplots()
ax3.set_xlabel('Volume in reactor [m^3]')
ax3.set_ylabel('Temperature [K]')
plt.grid()
plt.plot(plt_v,temperature,label="temperature")
plt.legend()
# plt.savefig(foldername+"Temperature.png")   
plt.show()


#%matplotlib inline #jupyter notebookでは必要
import numpy as np
from scipy.optimize import fsolve #熱バランス解くときに使った
import pandas #csvファイルの入出力系
import process_design_module2 as pdm #自作
from HYSYS_python_spreadsheets import Aspen_connection
import time

# 1.0 Data of the Aspen HYSYS file
File         = '20240508NRTL finished1.hsc'
Spreadsheets = ('SS_Reactor_outflow','SS_SeparateGas',"SS_Product")
Units        = ('E1', 'V1', 'VLV1', 'V2', 'E2', 
                'T1', 'E3', 'E4')

# 2.0 Perform connection
process_1      = Aspen_connection(File, Spreadsheets, Units)
# Test1.SSで次のやつができる　辞書型
#{'SS_Flash': <COMObject Item>, 'SS_turbine': <COMObject Item>, 'SS_Distillation': <COMObject Item>}

hy_ms = process_1.MaterialStreams
"""Item[0]=outflow"""
"""Item[2]=recycle"""




Product = process_1.SS['SS_Product']
Molar_product   = Product.Cell(1,0)             
y_meoh_product  = Product.Cell(1,1)

 
solver      = process_1.Solver

# 反応器の作成
reactor = pdm.Reactor(Pressure=30,) #標準では 論文と同じ条件で作成
P=3000
reactor.S = reactor.S * 2750.0  #4801個の反応器の並列分だけ断面積を変化
reactor.vmax = reactor.vmax * 2750.0 
a0 = 3 # [m s^-1] # 空塔速度
u0 = a0*reactor.S #[m^3/s]　#体積速度

# シミュレーションのステップ数の決定
dv = 0.01 #[m^3] # 微小体積
n_step = int(reactor.vmax/dv) #ステップ数
# 微小体積の流量の情報を保存

allowable_err = 0.00001
purge_ratio = 0.99

def objective_func(Ft):# 原料流量入れた際の反応器の組成出口誤差
    #シミュレーション
    F_co2 = Ft[0]
    F_h2 = Ft[1]
    print("F_co2[mol/s]:",F_co2," F_h2[mol/s]",F_h2)
    err = 1
    simulate_cnt = 1
    hy_ms.Item(0).Pressure.setValue(P)
    while (err > allowable_err) and (simulate_cnt < 10) :
        print(simulate_cnt)
        flows=list() 
        for n in range(n_step):
            flows.append(pdm.Flow())
        hy_tmp_frac = hy_ms.Item(2).ComponentMolarFraction.Values
        flows[0].Molar_flow["co"]   = purge_ratio*hy_ms.Item(2).MolarFlow.getValue("gmole/s")*hy_tmp_frac[3]
        flows[0].Molar_flow["co2"]  = F_co2
        flows[0].Molar_flow["h2"]   = F_h2
        flows[0].Molar_flow["h2o"]  = purge_ratio*hy_ms.Item(2).MolarFlow.getValue("gmole/s")*hy_tmp_frac[2]
        flows[0].Molar_flow["meoh"] = purge_ratio*hy_ms.Item(2).MolarFlow.getValue("gmole/s")*hy_tmp_frac[0]
        flows[0].Molar_flow["n2"]   = purge_ratio*hy_ms.Item(2).MolarFlow.getValue("gmole/s")*hy_tmp_frac[5]
        Temperature=np.full(n_step,485.0)
        for n in range(n_step-1):
            #反応速度式の計算
            xi1 = pdm.reaction_rate_1(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
            xi2 = pdm.reaction_rate_2(Temperature[n],reactor,flows[n])*reactor.catalyst*dv
            #モル流量の更新 F_out = F_in + (a*r1+b*r2)dv
            flows[n+1].Molar_flow["co"]   = flows[n].Molar_flow["co"]   + ( 0*xi1) + ( 1*xi2)
            flows[n+1].Molar_flow["co2"]  = flows[n].Molar_flow["co2"]  + (-1*xi1) + (-1*xi2)
            flows[n+1].Molar_flow["h2"]   = flows[n].Molar_flow["h2"]   + (-3*xi1) + (-1*xi2)
            flows[n+1].Molar_flow["h2o"]  = flows[n].Molar_flow["h2o"]  + ( 1*xi1) + ( 1*xi2)
            flows[n+1].Molar_flow["meoh"] = flows[n].Molar_flow["meoh"] + ( 1*xi1) + ( 0*xi2)
            flows[n+1].Molar_flow["n2"]   = flows[n].Molar_flow["n2"] #不活性
            #熱の変更
            def heat_balance(T):
                return pdm.Sum_FH(T,flows[n+1])- \
                        pdm.Sum_FH(Temperature[n],flows[n+1])+ \
                        pdm.Heat_of_reaction(Temperature[n],xi1,xi2)
            Temperature[n+1],*_ =fsolve(heat_balance,Temperature[n]) #heat_balanece(T)=0 となるTを求める
        #終了条件の確認
        hy_tmp_frac = hy_ms.Item(0).ComponentMolarFraction.Values
        err =(flows[-1].Molar_flow["meoh"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[0])**2 + \
            (flows[-1].Molar_flow["co2"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[1])**2 + \
            (flows[-1].Molar_flow["h2o"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[2])**2 + \
            (flows[-1].Molar_flow["co"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[3])**2 + \
            (flows[-1].Molar_flow["h2"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[4])**2 +\
            (flows[-1].Molar_flow["n2"] - hy_ms.Item(0).MolarFlow.getValue("gmole/s")*hy_tmp_frac[5])**2
        print("err =",err)
        data = list() #微小体積毎の温度・モル流量をcsvfileで保存
        for n in range(n_step): #各体積での処理
            data.append({"V[m^3]":dv*n,"T[K]":Temperature[n],"f_co[mol/s]":flows[n].Molar_flow["co"],"f_co2[mol/s]":flows[n].Molar_flow["co2"],"f_h2[mol/s]":flows[n].Molar_flow["h2"],"f_h2o[mol/s]":flows[n].Molar_flow["h2o"],"f_meoh[mol/s]":flows[n].Molar_flow["meoh"],"f_n2[mol/s]":flows[n].Molar_flow["n2"]})            
        data_csv = pandas.DataFrame(data)
        print(data_csv)
        simulate_cnt = simulate_cnt + 1
        solver.CanSolve = False
        hy_ms.Item(0).MolarFlow.setValue((flows[-1].Molar_flow["co"] + flows[-1].Molar_flow["co2"] + \
        flows[-1].Molar_flow["h2"] + flows[-1].Molar_flow["h2o"] +\
        flows[-1].Molar_flow["meoh"]+flows[-1].Molar_flow["n2"] )/(10**3)) #? mol/s
        hy_ms.Item(0).ComponentMolarFraction.Values = (flows[-1].y("meoh"),flows[-1].y("co2"),flows[-1].y("h2o"),flows[-1].y("co"),flows[-1].y("h2"),flows[-1].y("n2")) 
        hy_ms.Item(0).Temperature.setValue(Temperature[-1]-273.15) # セルシウス温度
        solver.CanSolve = True
        while solver.IsSolving == True: #解くまでまつ
            time.sleep(0.001) 
    print("finish.")

#試しにこれでやってみた
objective_func([2100.0,8000.0])
print("output")
# print(data_csv)

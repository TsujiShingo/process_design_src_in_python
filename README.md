# 更新履歴

- 最新版
    - process_design_moudule3.py：
    より汎用的なモジュール(反応の追加などをしやすくしました)
    - process_design_simulate3.py：
    CUIベースで等温反応器、非等温断熱反応器を扱います。
    現時点ではHysysとつなげません。
    moudule3.pyが必要です。
    - (非公開、研究室内公開) main.py 
    上記2つのモジュールと
    https://github.com/edgarsmdn/Aspen_HYSYS_Python.git
    にあるHYSYS_python_spreadsheets.pyを用いて、
    CUIベースで等温反応器、非等温断熱反応器を扱います。
    Hysysと接続しなくても一部の機能は利用可能です。
    

- そのほか
    - process_design_moudule.py：
    一番最初に作ったモジュールです
    - process_design_moudule2.py：
        排ガスに含まれる不活性物質として窒素を入れたものです。
    - processdesign_2024_kucpe_12_course_N2_ver1.py：
        排ガスに含まれる不活性物質として窒素を入れたものです。
        module2.pyが必要です
    - processdesign_2024_kucpe_12_course_version1.py：
        一番最初に試作で作ったシミュレーションです。
        ここでのmodule.pyは必要ありません。
    - processdesign_2024_kucpe_12_course_version2.py：
        version1.pyを更新したものです。
        ここでのmodule.pyは必要ありません。
    - processdesign_2024_kucpe_12_course_version3_module.py
        version2.pyのコードが長くなったので、フローや反応器、熱容量の計算などをmodule.pyに移しました。
        module.pyが必要です。
    - processdesign_2024_kucpe_12_course_version4.py
        Aspen Hysysとつなげて原料流量一定の際の収束計算をするものです。
        module.pyの他に、
        https://github.com/edgarsmdn/Aspen_HYSYS_Python.git
        にあるHYSYS_python_spreadsheets.pyが必要です。

    - processdesign_2024_kucpe_12_course_version5.py
        Aspen Hysysとつなげて目標生産物が得られるような原料流量を計算します。
        module.pyの他に、
        https://github.com/edgarsmdn/Aspen_HYSYS_Python.git
        にあるHYSYS_python_spreadsheets.pyが必要です。


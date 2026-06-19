# XELYM
XELYM (XlinkEd poLYmer gel Modeling tools)のソースコード（開発段階中）

基本的にjson形式のファイルを読み込ます形式
コンパイル後，例えば以下のようなjsonのインプットファイルを用意した上で，
{
"mode": "crosslinking-radical",
"trajfile": "../02_nvt/nvt.xtc",
"itpfile": [
 "../sys/nipam_bis.itp",
 "../data_set/itp/imit.itp"
],
"outitp": "nipam_bis_temp.itp",
"chargefile": [
 "../data_set/charge/nipam/charge_nipam.dat",
 "../data_set/charge/bis/charge_0-1-0-0.dat",
 "../data_set/charge/bis/charge_0-0-0-1.dat",
 "../data_set/charge/bis/charge_0-1-0-1.dat",
 "../data_set/charge/bis/charge_1-1-0-0.dat",
 "../data_set/charge/bis/charge_0-0-1-1.dat",
 "../data_set/charge/bis/charge_0-1-1-1.dat",
 "../data_set/charge/bis/charge_1-1-0-1.dat",
 "../data_set/charge/bis/charge_1-1-1-1.dat"],
"selpol": ["C1_Ni", "C2_Ni"],
"selcross": ["C2_BIS", "C1_BIS", "C6_BIS", "C7_BIS"],
"startlist": [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600],
"Np": 600,
"Nc": 60,
"rc": 4.0
}
./xelym.x run.json
と実行することで，自動的にモデリングすることは可能

まだ開発段階なので，マニュアルなどは追って作成する予定である．

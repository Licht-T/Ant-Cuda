# アリさん採餌シミュレーションマニュアル

## 実行方法

1. `freeglut`とかをインストール（LinuxだとたいていCUDA入れたときに入る）．
2. `batch.py`の`dists`リストにチェックしたい巣と餌の距離を，`angles`リストにチェックしたい巣と2つの餌の角度を追加する．
2. `python ./batch.py`するとコンパイルが実行され計算が始まる．データはカレントディレクトリに新たに作成されるフォルダにセーブされる．

## より詳しい説明

### 定数に関して

（おおむね）すべての定数は`Constants.h`に定義されたマクロで表される．
（入ってないのは，Critical Angleの判定定数とヒルの式のシフトパラメタぐらいか？）

* `MACRO_NMAX`: アリの数
* `MACRO_MAX`: セルの縦横サイズ
* `MACRO_FOODSOURCE`: 餌の初期量
* `MACRO_NUM_FOODS`: 餌の数
* `MACRO_FOOD_DIST`: 巣と餌の距離
* `MACRO_FOOD_ANGLE`: 巣と2つの餌のなす角度
* `MACRO_NEST_X`: 巣のx座標（上からMACRO_NEST_X列目のセル）
* `MACRO_NEST_Y`: 巣のy座標（左からMACRO_NEST_Y行目のセル）
* `MACRO_REC`: 餌のRecovery Ratio
* `MACRO_EVAPOLATION_CONST`: フェロモンの蒸発率
* `MACRO_MAX_SEARCH_TIME`: アリの最長餌探索時間（この時間を過ぎると強制的に帰巣モードに）
* `MACRO_MAX_STEP`: シミュレーションの試行回数（この回数シミュレーションをやって採餌効率のアンサンブル平均をだす）
* `MACRO_MAX_TIME`: 1シミュレーションの実行時間の最大値（tがこの時間に到達したら1シミュレーションが終了）
* `MACRO_OFFSET_TIME`: 採餌効率をだすときに無視する初期状態の時間数
* `MACRO_EMI`: フェロモンの放出量
* `MACRO_ENEST`: 紆余曲折を経て残ってしまっている定数だが，ただ1倍しているだけなので問題ない（）
* `MACRO_HIL_CONST`: ヒルの式の定数alpha
* `MACRO_UNIT`: 採餌一回あたり減る餌の量
* `MACRO_DIFFE`: 拡散定数

### グリッドサーチの方法

CUDAデバイスとホスト・コンピュータの間のデータ齟齬を避けるため，全ての実験時定数をマクロで与えている．
そのため，グリッドサーチしたいパラメタの組み合わせをビルド時にマクロコマンド（`-D...`コマンド）として与え，
それらの組み合わせごとにバイナリを作ることでグリッドサーチを実現する．

1. グリッドサーチしたいパラメタを`Constants.h`からコメントアウト．
2. `batch.py`をうまく書き換えて`make`時にグリッドサーチをしたいパラメタを引数として与えてやるように設定．
3. `makefile`をうまく書き換えて`make`の引数として渡ってきたパラメタをマクロ定数としてビルドするように設定．
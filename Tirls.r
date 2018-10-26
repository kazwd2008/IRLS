#####################################################################
#  繰返し再加重最小二乗法(IRLS)によるM推定 ウエイトはTukeyのbiweight#
#-------------------------------------------------------------------#
#       Ver. 0    2010/06/14  オリジナル                            #
#       Ver. 1    2010/12/01  乗率対応版                            #
#       Ver. 1.1  2011/03/07  無限ループ回避                        #
#       Ver. 2.0  2012/07/19  コンスタントデフォルト値修正          #
#       Ver. 3.0  2012/12/18  平均絶対偏差版と中央絶対偏差版作成    #
#                             Tirls.ave       Tirls.med             #
#       Ver. 3.1  2012/12/19  収束条件パラメータ追加                #
#       Ver. 3.2  2013/04/25  切片のない回帰版を作成
#                             Tirls.avex     Tirls.medxを作成  					#
#       Ver. 3.3  2014/08/05  WLSの場合残差平均0ではないので        #
#                               AADの計算を修正 <= 間違い 09/05		#
#       Ver. 3.4  2014/08/25  初期値計算はOLS                       #
#       Ver. 4.0  2014/09/05  残差から平均・中央値を引かない        #
#       Ver. 4.1  2018/10/26  公開用にコメント整備
#####################################################################
#  収録:  Tukeyのbiweightを使用したIRLSアルゴリズムの関数
#    Tirls.ave     尺度AAD, 回帰モデル用
#    Tirls.med     尺度MAD, 回帰モデル用
#    Tirls.avex    尺度AAD, 切片のない回帰モデル
#    Tirls.medx   尺度MAD,  切片のない回帰モデル
#    RTirls.ave    尺度AAD, 比率モデル
#####################################################################
#  関数パラメータ
#   y1      目的変数（必須）                               
#   x1      説明変数（必須）
#   rt      乗率。指定がない場合デフォルトは1
#   c1      biweight関数用コンスタント。デフォルトはgeneralにパフォーマ
#            ンスが良い最大値(aveで8)。 
#   c=4     とてもロバスト 
#   c=8     ちょっとロバスト
#   dat     x1, y1が含まれるデータフレーム。指定がなければ使わない
#   rp.max  ループ回数の上限。特に指定しなければ150回  
###################################################################
#  関数戻り値
#   TK      最終結果
#   wt      ウエイト
#   rp      ループ回数
#   s1      平均絶対残差(AAD)または中央絶対偏差(MAD)の変遷データ
###################################################################

Tirls.ave <- function(y1, x1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # データフレーム指定があるときのみ
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1の変遷を保存
  R0 <- lm(y1~x1, weights=rt) 
  #R0 <- lm(y1~x1)       				# 初期値は普通のOLS
  
  Tk2 <- R0
  rp1 <- 1								# ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))       # 平均絶対残差(MAD)

  #### ウエイト算出
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ループカウンタ
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# 平均絶対残差
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.med <- function(y1, x1, rt=rep(1, length(y1)), c1=10.03, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # データフレーム指定があるときのみ
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1, weights=rt) 
  #R0 <- lm(y1~x1)       				# 初期値は普通のOLS
  
  Tk2 <- R0
  rp1 <- 1								# ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mad(Tk2$residuals) 			# 中央絶対残差(MAD) σ準拠に補正済み

  #### ウエイト算出
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.18
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ループカウンタ
        s1 <- s1.cg[rp1] <- mad(Tk2$residuals)	# 中央絶対残差(MAD)
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.avex <- function(y1, x1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # データフレーム指定があるときのみ
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1の変遷を保存
  R0 <- lm(y1~x1-1, weights=rt)  
  #R0 <- lm(y1~x1-1)       				# 初期値は普通のOLS
  
  Tk2 <- R0
  rp1 <- 1								# ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))   # 平均絶対残差(AAD)

  #### ウエイト算出
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1-1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ループカウンタ
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# 平均絶対残差
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.medx <- function(y1, x1, rt=rep(1, length(y1)), c1=10.3, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # データフレーム指定があるときのみ
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1-1, weights=rt)  
  #R0 <- lm(y1~x1-1)       				# 初期値は普通のOLS
  
  Tk2 <- R0
  rp1 <- 1								# ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mad(Tk2$residuals) # 中央絶対残差(MAD) σ準拠に補正済み

  #### ウエイト算出
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.18
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1-1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ループカウンタ
        s1 <- s1.cg[rp1] <- mad(Tk2$residuals)		# 中央絶対残差(MAD)
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
RTirls.ave <- function(y1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # データフレーム指定があるときのみ
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1の変遷を保存
  R0 <- lm(y1~1, weights=rt)  
  #R0 <- lm(y1~1)       				# 初期値は普通のOLS
  
  Tk2 <- R0
  rp1 <- 1								# ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))   # 平均絶対残差(AAD)

  #### ウエイト算出
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ループカウンタ
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# 平均絶対残差
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------


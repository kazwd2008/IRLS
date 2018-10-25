#####################################################################
#  繰返し再加重最小二乗法(IRLS)によるM推定  Huberウエイト版         #
#-------------------------------------------------------------------#
#       Ver. 0    2010/08/05  オリジナル                            #
#       Ver. 1    2010/12/01  乗率対応版                            #
#       Ver. 1.1  2011/03/07  大量ループ回避                        #
#       Ver. 2.0  2012/07/19  コンスタントc1修正                    #
#       Ver. 3.0  2012/12/18  平均絶対偏差版と中央絶対偏差版作成    #
#                             Hirls.ave       Hirls.med             #
#       Ver. 3.1  2012/12/19  収束条件パラメータ追加                #
#       Ver. 3.2  2013/04/25  切片のない回帰用のHirls.avexと        #
#                               Hirls.madxを作成                    #
#       Ver. 3.3  2014/08/05  WLSの場合残差平均0ではないので        #
#                               AADの計算を修正  					#
#####################################################################
#  関数Hirls パラメータ
#   y1      単回帰の目的変数（必須）
#   x1      単回帰の説明変数（必須）
#   rt      乗率。指定がない場合デフォルトは1
#   c1      ウエイト関数用コンスタント。デフォルトは1.686(正規分布想定)。
#       1.44でTukeyのbiweight c=4、 2.16でc=6、2.88でTukeyのc=8相当。
#   dat     x1, y1が含まれるデータフレーム。指定がなければ使わない
#   rp.max  ループ回数の上限。特に指定しなければ50回  
###################################################################
#  関数Hirls 戻り値
#   HB      最終結果
#   wt      ウエイト
#   rp      ループ回数
#   s1      平均絶対残差(AAD)または中央絶対偏差(MAD)の変遷データ
###################################################################
#	Hirls.avexとHirls.madxは、切片のない回帰を行う。
###################################################################

Hirls.ave <- function(y1, x1, rt=rep(1, length(y1)), c1=1.15, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # データフレーム指定があるときのみ 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1, weights=rt)						            # 初期値はウエイトをつけないOLS
  
  Hb2 <- R0
  rp1 <- 1                              # ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))        # 平均絶対残差

  #### ウエイト算出
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ループカウンタ
      s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))	# 平均絶対残差
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Hirls.med <- function(y1, x1, rt=rep(1, length(y1)), c1=1.44, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # データフレーム指定があるときのみ 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1, weights=rt)						            # 初期値はウエイトをつけないOLS
  
  Hb2 <- R0
  rp1 <- 1                              # ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mad(Hb2$residuals)        # 中央絶対残差  σ準拠に補正済み

  #### ウエイト算出
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ループカウンタ
      s1 <- s1.cg[rp1] <- mad(Hb2$residuals)	# 中央絶対残差
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------

Hirls.avex <- function(y1, x1, rt=rep(1, length(y1)), c1=1.15, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # データフレーム指定があるときのみ 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1-1, weights=rt)						            # 初期値はウエイトをつけないOLS
  
  Hb2 <- R0
  rp1 <- 1                              # ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))        # 平均絶対残差

  #### ウエイト算出
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1-1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ループカウンタ
      s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))	# 平均絶対残差
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Hirls.medx <- function(y1, x1, rt=rep(1, length(y1)), c1=1.44, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # データフレーム指定があるときのみ 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MADの変遷を保存
  R0 <- lm(y1~x1-1, weights=rt)						            # 初期値はウエイトをつけないOLS
  
  Hb2 <- R0
  rp1 <- 1                              # ループ回数
  s0 <- 0			# 最低１回はループするように差が0.01より大きくなるよう値をセット
  s1 <- s1.cg[rp1] <- mad(Hb2$residuals)        # 中央絶対残差  σ準拠に補正済み

  #### ウエイト算出
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### 繰返し計算
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1-1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ループカウンタ
      s1 <- s1.cg[rp1] <- mad(Hb2$residuals)	# 中央絶対残差
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------

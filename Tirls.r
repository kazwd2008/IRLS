#####################################################################
#  �J�Ԃ��ĉ��d�ŏ����@(IRLS)�ɂ��M���� �E�G�C�g��Tukey��biweight#
#-------------------------------------------------------------------#
#       Ver. 0    2010/06/14  �I���W�i��                            #
#       Ver. 1    2010/12/01  �旦�Ή���                            #
#       Ver. 1.1  2011/03/07  �������[�v���                        #
#       Ver. 2.0  2012/07/19  �R���X�^���g�f�t�H���g�l�C��          #
#       Ver. 3.0  2012/12/18  ���ϐ�Ε΍��łƒ�����Ε΍��ō쐬    #
#                             Tirls.ave       Tirls.med             #
#       Ver. 3.1  2012/12/19  ���������p�����[�^�ǉ�                #
#       Ver. 3.2  2013/04/25  �ؕЂ̂Ȃ���A�ł��쐬
#                             Tirls.avex     Tirls.medx���쐬  					#
#       Ver. 3.3  2014/08/05  WLS�̏ꍇ�c������0�ł͂Ȃ��̂�        #
#                               AAD�̌v�Z���C�� <= �ԈႢ 09/05		#
#       Ver. 3.4  2014/08/25  �����l�v�Z��OLS                       #
#       Ver. 4.0  2014/09/05  �c�����畽�ρE�����l�������Ȃ�        #
#       Ver. 4.1  2018/10/26  ���J�p�ɃR�����g����
#####################################################################
#  ���^:  Tukey��biweight���g�p����IRLS�A���S���Y���̊֐�
#    Tirls.ave     �ړxAAD, ��A���f���p
#    Tirls.med     �ړxMAD, ��A���f���p
#    Tirls.avex    �ړxAAD, �ؕЂ̂Ȃ���A���f��
#    Tirls.medx   �ړxMAD,  �ؕЂ̂Ȃ���A���f��
#    RTirls.ave    �ړxAAD, �䗦���f��
#####################################################################
#  �֐��p�����[�^
#   y1      �ړI�ϐ��i�K�{�j                               
#   x1      �����ϐ��i�K�{�j
#   rt      �旦�B�w�肪�Ȃ��ꍇ�f�t�H���g��1
#   c1      biweight�֐��p�R���X�^���g�B�f�t�H���g��general�Ƀp�t�H�[�}
#            ���X���ǂ��ő�l(ave��8)�B 
#   c=4     �ƂĂ����o�X�g 
#   c=8     ������ƃ��o�X�g
#   dat     x1, y1���܂܂��f�[�^�t���[���B�w�肪�Ȃ���Ύg��Ȃ�
#   rp.max  ���[�v�񐔂̏���B���Ɏw�肵�Ȃ����150��  
###################################################################
#  �֐��߂�l
#   TK      �ŏI����
#   wt      �E�G�C�g
#   rp      ���[�v��
#   s1      ���ϐ�Ύc��(AAD)�܂��͒�����Ε΍�(MAD)�̕ϑJ�f�[�^
###################################################################

Tirls.ave <- function(y1, x1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # �f�[�^�t���[���w�肪����Ƃ��̂�
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1�̕ϑJ��ۑ�
  R0 <- lm(y1~x1, weights=rt) 
  #R0 <- lm(y1~x1)       				# �����l�͕��ʂ�OLS
  
  Tk2 <- R0
  rp1 <- 1								# ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))       # ���ϐ�Ύc��(MAD)

  #### �E�G�C�g�Z�o
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ���[�v�J�E���^
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# ���ϐ�Ύc��
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.med <- function(y1, x1, rt=rep(1, length(y1)), c1=10.03, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # �f�[�^�t���[���w�肪����Ƃ��̂�
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1, weights=rt) 
  #R0 <- lm(y1~x1)       				# �����l�͕��ʂ�OLS
  
  Tk2 <- R0
  rp1 <- 1								# ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mad(Tk2$residuals) 			# ������Ύc��(MAD) �Џ����ɕ␳�ς�

  #### �E�G�C�g�Z�o
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.18
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ���[�v�J�E���^
        s1 <- s1.cg[rp1] <- mad(Tk2$residuals)	# ������Ύc��(MAD)
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.avex <- function(y1, x1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # �f�[�^�t���[���w�肪����Ƃ��̂�
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1�̕ϑJ��ۑ�
  R0 <- lm(y1~x1-1, weights=rt)  
  #R0 <- lm(y1~x1-1)       				# �����l�͕��ʂ�OLS
  
  Tk2 <- R0
  rp1 <- 1								# ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))   # ���ϐ�Ύc��(AAD)

  #### �E�G�C�g�Z�o
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1-1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ���[�v�J�E���^
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# ���ϐ�Ύc��
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Tirls.medx <- function(y1, x1, rt=rep(1, length(y1)), c1=10.3, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # �f�[�^�t���[���w�肪����Ƃ��̂�
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1-1, weights=rt)  
  #R0 <- lm(y1~x1-1)       				# �����l�͕��ʂ�OLS
  
  Tk2 <- R0
  rp1 <- 1								# ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mad(Tk2$residuals) # ������Ύc��(MAD) �Џ����ɕ␳�ς�

  #### �E�G�C�g�Z�o
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.18
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~x1-1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ���[�v�J�E���^
        s1 <- s1.cg[rp1] <- mad(Tk2$residuals)		# ������Ύc��(MAD)
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
RTirls.ave <- function(y1, rt=rep(1, length(y1)), c1=8, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))	        # �f�[�^�t���[���w�肪����Ƃ��̂�
  s1.cg <- rep(0, rp.max)               # 2012.12.18 s1�̕ϑJ��ۑ�
  R0 <- lm(y1~1, weights=rt)  
  #R0 <- lm(y1~1)       				# �����l�͕��ʂ�OLS
  
  Tk2 <- R0
  rp1 <- 1								# ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))   # ���ϐ�Ύc��(AAD)

  #### �E�G�C�g�Z�o
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- lm(y1~1, weights=w1*rt) 
        rp1 <- rp1 + 1			# ���[�v�J�E���^
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# ���ϐ�Ύc��
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------


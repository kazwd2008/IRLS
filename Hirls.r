#####################################################################
#  �J�Ԃ��ĉ��d�ŏ����@(IRLS)�ɂ��M����  Huber�E�G�C�g��         #
#-------------------------------------------------------------------#
#       Ver. 0    2010/08/05  �I���W�i��                            #
#       Ver. 1    2010/12/01  �旦�Ή���                            #
#       Ver. 1.1  2011/03/07  ��ʃ��[�v���                        #
#       Ver. 2.0  2012/07/19  �R���X�^���gc1�C��                    #
#       Ver. 3.0  2012/12/18  ���ϐ�Ε΍��łƒ�����Ε΍��ō쐬    #
#                             Hirls.ave       Hirls.med             #
#       Ver. 3.1  2012/12/19  ���������p�����[�^�ǉ�                #
#       Ver. 3.2  2013/04/25  �ؕЂ̂Ȃ���A�p��Hirls.avex��        #
#                               Hirls.madx���쐬                    #
#       Ver. 3.3  2014/08/05  WLS�̏ꍇ�c������0�ł͂Ȃ��̂�        #
#                               AAD�̌v�Z���C��  					#
#####################################################################
#  �֐�Hirls �p�����[�^
#   y1      �P��A�̖ړI�ϐ��i�K�{�j
#   x1      �P��A�̐����ϐ��i�K�{�j
#   rt      �旦�B�w�肪�Ȃ��ꍇ�f�t�H���g��1
#   c1      �E�G�C�g�֐��p�R���X�^���g�B�f�t�H���g��1.686(���K���z�z��)�B
#       1.44��Tukey��biweight c=4�A 2.16��c=6�A2.88��Tukey��c=8�����B
#   dat     x1, y1���܂܂��f�[�^�t���[���B�w�肪�Ȃ���Ύg��Ȃ�
#   rp.max  ���[�v�񐔂̏���B���Ɏw�肵�Ȃ����50��  
###################################################################
#  �֐�Hirls �߂�l
#   HB      �ŏI����
#   wt      �E�G�C�g
#   rp      ���[�v��
#   s1      ���ϐ�Ύc��(AAD)�܂��͒�����Ε΍�(MAD)�̕ϑJ�f�[�^
###################################################################
#	Hirls.avex��Hirls.madx�́A�ؕЂ̂Ȃ���A���s���B
###################################################################

Hirls.ave <- function(y1, x1, rt=rep(1, length(y1)), c1=1.15, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # �f�[�^�t���[���w�肪����Ƃ��̂� 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1, weights=rt)						            # �����l�̓E�G�C�g�����Ȃ�OLS
  
  Hb2 <- R0
  rp1 <- 1                              # ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))        # ���ϐ�Ύc��

  #### �E�G�C�g�Z�o
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ���[�v�J�E���^
      s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))	# ���ϐ�Ύc��
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Hirls.med <- function(y1, x1, rt=rep(1, length(y1)), c1=1.44, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # �f�[�^�t���[���w�肪����Ƃ��̂� 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1, weights=rt)						            # �����l�̓E�G�C�g�����Ȃ�OLS
  
  Hb2 <- R0
  rp1 <- 1                              # ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mad(Hb2$residuals)        # ������Ύc��  �Џ����ɕ␳�ς�

  #### �E�G�C�g�Z�o
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ���[�v�J�E���^
      s1 <- s1.cg[rp1] <- mad(Hb2$residuals)	# ������Ύc��
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------

Hirls.avex <- function(y1, x1, rt=rep(1, length(y1)), c1=1.15, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # �f�[�^�t���[���w�肪����Ƃ��̂� 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1-1, weights=rt)						            # �����l�̓E�G�C�g�����Ȃ�OLS
  
  Hb2 <- R0
  rp1 <- 1                              # ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))        # ���ϐ�Ύc��

  #### �E�G�C�g�Z�o
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1-1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ���[�v�J�E���^
      s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals-mean(Hb2$residuals)))	# ���ϐ�Ύc��
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-----------------------------------------------------------------------------
Hirls.medx <- function(y1, x1, rt=rep(1, length(y1)), c1=1.44, dat="", rp.max=150, cg.rt=0.01) {

  if (dat!="") attach(get(dat))         # �f�[�^�t���[���w�肪����Ƃ��̂� 
  s1.cg <- rep(0, rp.max)               # 2011.03.07 MAD�̕ϑJ��ۑ�
  R0 <- lm(y1~x1-1, weights=rt)						            # �����l�̓E�G�C�g�����Ȃ�OLS
  
  Hb2 <- R0
  rp1 <- 1                              # ���[�v��
  s0 <- 0			# �Œ�P��̓��[�v����悤�ɍ���0.01���傫���Ȃ�悤�l���Z�b�g
  s1 <- s1.cg[rp1] <- mad(Hb2$residuals)        # ������Ύc��  �Џ����ɕ␳�ς�

  #### �E�G�C�g�Z�o
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### �J�Ԃ��v�Z
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- lm(y1~x1-1, weights=w1*rt) 
      rp1 <- rp1 + 1			# ���[�v�J�E���^
      s1 <- s1.cg[rp1] <- mad(Hb2$residuals)	# ������Ύc��
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-----------------------------------------------------------------------------

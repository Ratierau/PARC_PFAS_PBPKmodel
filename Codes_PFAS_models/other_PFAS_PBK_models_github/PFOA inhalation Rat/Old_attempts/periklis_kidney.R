
CKBf <- solution$CKBf
CKFf <-  solution$CKFf
CKTf <- solution$CKTf

in_blood <- CKBf*QBK/500 
out_lymph <- - CKFf*QBK/500
dif_blood <- PeffK*AK*(CKBf-CKFf)
dif_tis<- - kKFKT*(CKFf-CKTf)
out_ota1<-    -(VmK_Oat1*CKFf/KmK_Oat1+CKFf)*VKF
out_oat2 <- - (VmK_Oat3*CKFf/KmK_Oat3+CKFf)*VKF
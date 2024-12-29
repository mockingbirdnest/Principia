PUBLIC identity@@8
PUBLIC sqrtps_xmm0_xmm0@@8
PUBLIC mulsd_xmm0_xmm0@@8
PUBLIC mulsd_xmm0_xmm0_4x@@8
PUBLIC sqrtsd_xmm0_xmm0@@8

.CODE
identity@@8 PROC
  ret
identity@@8 ENDP

sqrtps_xmm0_xmm0@@8 PROC
  sqrtps xmm0, xmm0
  ret
sqrtps_xmm0_xmm0@@8 ENDP

mulsd_xmm0_xmm0@@8 PROC
  mulsd xmm0, xmm0
  ret
mulsd_xmm0_xmm0@@8 ENDP

mulsd_xmm0_xmm0_4x@@8 PROC
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  ret
mulsd_xmm0_xmm0_4x@@8 ENDP

sqrtsd_xmm0_xmm0@@8 PROC
  sqrtsd xmm0, xmm0
  ret
sqrtsd_xmm0_xmm0@@8 ENDP
END

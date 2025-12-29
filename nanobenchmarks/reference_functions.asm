PUBLIC identity@@8
PUBLIC sqrtps_xmm0_xmm0@@8
PUBLIC mulsd_xmm0_xmm0@@8
PUBLIC mulsd_xmm0_xmm0_4x@@8
PUBLIC sqrtsd_xmm0_xmm0@@8
PUBLIC fill3_from_xmm0
PUBLIC fill6_from_xmm0
PUBLIC fill_from_ymm0
PUBLIC fill_from_ymm0_ymm1

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

fill3_from_xmm0 PROC
  vxorpd      xmm1, xmm1, xmm1
  vunpcklpd   xmm2, xmm0, xmm1
  vmovddup    xmm1, xmm2
  vinsertf128 ymm2, ymm1, xmm2, 1
  vmovupd     YMMWORD PTR [rcx], ymm2
  mov         rax, rcx
  vzeroupper
  ret
fill3_from_xmm0 ENDP

fill_from_ymm0 PROC
  vmovupd     YMMWORD PTR [rcx], ymm0
  mov         rax, rcx
  vzeroupper
  ret
fill_from_ymm0 ENDP

fill6_from_xmm0 PROC
  vxorpd      xmm1, xmm1, xmm1
  vunpcklpd   xmm2, xmm0, xmm1
  vmovddup    xmm1, xmm2
  vinsertf128 ymm2, ymm1, xmm2, 1
  vmovupd     YMMWORD PTR [rcx], ymm2
  vmovupd     YMMWORD PTR [rcx+32], ymm2
  mov         rax, rcx
  vzeroupper
  ret
fill6_from_xmm0 ENDP

fill_from_ymm0_ymm1 PROC
  vmovupd     YMMWORD PTR [rcx], ymm0
  vmovupd     YMMWORD PTR [rcx+32], ymm1
  mov         rax, rcx
  vzeroupper
  ret
fill_from_ymm0_ymm1 ENDP

END

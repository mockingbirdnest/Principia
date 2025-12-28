PUBLIC foofoofoo
PUBLIC googoogoo

.CODE
foofoofoo PROC
  vxorpd      xmm1, xmm1, xmm1
  vunpcklpd   xmm2, xmm0, xmm1
  vmovddup    xmm1, xmm2
  vinsertf128 ymm2, ymm1, xmm2, 1
  vmovupd     YMMWORD PTR [rcx], ymm2
  mov         rax, rcx
  vzeroupper
  ret
foofoofoo ENDP

googoogoo PROC
  vmovupd     YMMWORD PTR [rcx], ymm0
  mov         rax, rcx
  vzeroupper
  ret
googoogoo ENDP
END
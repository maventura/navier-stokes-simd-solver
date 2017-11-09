[BITS 64]

extern printf
global vvp_asm


section .data
copy_float_value_mask: DB 0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3
eltest: DB 0x05
wt: DD 0.8
uno: DD 1.0, 1.0, 1.0, 1.0
seis: DD 6.0, 6.0, 6.0, 6.0
dos: DD 2.0, 2.0, 2.0, 2.0

section .text
%define q_ xmm11
%define h_ xmm12
%define wt_ xmm13
%define delta_ xmm14
%define r_ xmm15
%define offset_i_ rdx
%define offset_j_ rcx
%define offset_k_ 1
%define mat_arr_ rdi
%define pos_ rsi



%define offset_U2_ 0
%define offset_V2_ 8
%define offset_W2_ 16
%define offset_omx2_ 24
%define offset_omy2_ 32
%define offset_omz2_ 40
%define offset_psix2_ 48
%define offset_psiy2_ 56
%define offset_psiz2_ 64
%define offset_omx1_ 72
%define offset_omy1_ 80
%define offset_omz1_ 88
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

vvp_asm:	

	push rbp
	mov rbp, rsp
	push r15  ;estos son los 5 registros a preservar... si es que los pisamos
	push r14
	push r13
	push r12
	push rbx
	sub rsp, 8  ;la pila debe estar alineada a 16 bits


  ;int vvp_asm  (
    ;float* data,     rdi
    ;int    pos       rsi
    ;float  r		  xmm0
    ;float  dx		  xmm1
    ;float  q		  xmm2
    ;int 	offset i  rdx
    ;int    offset j  rcx
  ;)

	xorps r_, r_
	movdqu r_, xmm0
	movdqu xmm0, [copy_float_value_mask]
	pshufb r_, xmm0

	xorps h_, h_
	movdqu h_, xmm1
	movdqu xmm1, [copy_float_value_mask]
	pshufb h_, xmm1

	xorps q_, q_
	movdqu q_, xmm2
	movdqu xmm2, [copy_float_value_mask]
	pshufb q_, xmm2

	xorps wt_, wt_
	movdqu wt_, [wt]
	movdqu xmm0, [copy_float_value_mask]
	pshufb wt_, xmm0

	; --------------------Calculo de eje x----------------------------------

	call calcular_v 

	;calculo de delta: (1 + q * [U2.at(i - 1, j, k) - U2.at(i + 1, j, k)] + 6 * r)
	mov rax, pos_
	sub rax, offset_i_ ; aca tenemos pos para i-1, j, k
	mov r15, [mat_arr_ + offset_U2_]
	movdqu xmm0, [r15 + rax * 4]
	mov rax, pos_
	add rax, offset_i_ ; aca tenemos pos para i+1, j, k
	movdqu xmm1, [r15 + rax * 4]
	subps xmm0, xmm1
	mulps xmm0, q_
	movdqu xmm2, [uno]
	addps xmm0, xmm2
	movdqu xmm2, [seis]
	mulps xmm2, r_
	addps xmm0, xmm2 
	movdqu delta_, xmm0

	;from now on xmm0 acumulates the result
	xorps xmm0, xmm0
	
	;p1 = (-U2.at(i + 1, j, k) + r)
	mov r15, [mat_arr_ + offset_U2_]
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1 ;TODO: XMM1 tiene que ser cero para hacer U2 negativo, era para eso no?
	subps xmm1, xmm2
	addps xmm1, r_
	;p1 * omx2.at(i + 1, j, k)
	mov r14, [mat_arr_ + offset_omx2_]
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1 

	;p3 = (U2.at(i - 1, j, k) + r)
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p3 * omx2.at(i - 1, j, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p2 = (-V2.at(i, j + 1, k) + r)
	mov r15, [mat_arr_ + offset_V2_]
	mov rax, pos_
	add rax, offset_j_
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1 ; TODO: IDEM que antes, xmm1 debe ser cero
	subps xmm1, xmm2
	addps xmm1, r_
	;p2 * omx2.at(i, j + 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p4 = (V2.at(i, j - 1, k) + r)
	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p4 * omx2.at(i, j - 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p5 = (W2.at(i, j, k + 1) + r)
	mov r15, [mat_arr_ + offset_W2_]
	mov rax, pos_
	add rax, offset_k_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p5 * omx2.at(i, j, k + 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p6 = (W2.at(i, j, k - 1) + r)
	mov rax, pos_
	sub rax, offset_k_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p6 * omx2.at(i, j, k - 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;(1.0 / q)*omx1.at(i, j, k)
	mov r15, [mat_arr_ + offset_omx1_]
	movdqu xmm1, [r15 + pos_*4]
	divps xmm1, q_
	addps xmm0, xmm1

	;ax1 = omy2.at(i, j, k) * U2.at(i, j + 1, k); 
	mov r15, [mat_arr_ + offset_omy2_]
	movdqu xmm1, [r15 + pos_*4]		;xmm1 referencia a omy2(i,j,k)
	mov r14, [mat_arr_ + offset_U2_]
	mov rax, pos_
	add rax, offset_j_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

    ;ax2 = -omy2.at(i, j, k) * U2.at(i, j - 1, k);
	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;ax3 = omz2.at(i, j, k) * U2.at(i, j, k + 1);
    mov r15, [mat_arr_ + offset_omz2_]
	movdqu xmm1, [r15 + pos_*4]		;xmm1 referencia a omz2(i,j,k)
	mov rax, pos_
	add rax, offset_k_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

    ;ax4 = -omz2.at(i, j, k) * U2.at(i, j, k - 1);
	mov rax, pos_
	sub rax, offset_k_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;* wt * (q / delta)
	mulps xmm0, q_
	divps xmm0, delta_
	mulps xmm0, wt_

	;(1.0 - wt)*omx2(i,j,k)_old 
    mov r15, [mat_arr_ + offset_omx2_]
	movdqu xmm1, [r15 + pos_*4]
	movdqu xmm2, [uno]
	subps xmm2, wt_
	mulps xmm1, xmm2
	addps xmm0, xmm1

    mov r15, [mat_arr_ + offset_omx2_]
	mov r14, pos_
	shl r14, 2
	add r15, r14
	movdqu [r15], xmm0 ;hasta aca calculamos omx2(i,j,k)


;------------------------Eliptica eje x----------------------------
	mov r15, [mat_arr_ + offset_psix2_]
	mov r14, pos_
	add r14, offset_i_
	movdqu xmm0, [r15 + r14*4] ;psix2.at(i + 1, j, k)
 
	mov r14, pos_
	add r14, offset_j_
	movdqu xmm3,  [r15 + r14*4]
	addps xmm0, xmm3;+ psix2.at(i, j + 1, k)

	mov r14, pos_
	sub r14, offset_i_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3 ;+ psix2.at(i - 1, j, k)
	;TODO: addps breaks if not aligned to 16.
	;decide if either replacing with two step addps or aligning the data
	;is the best option.

	mov r14, pos_
	sub r14, offset_j_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3 ;+ psix2.at(i, j - 1, k)

	mov r14, pos_
	add r14, offset_k_
	movdqu xmm3, [r15 + r14*4] 
	addps xmm0, xmm3;+ psix2.at(i, j, k + 1)

	mov r14, pos_
	sub r14, offset_k_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0,  xmm3;+ psix2.at(i, j, k - 1)

	mov r15, [mat_arr_ + offset_omx2_]
	movdqu xmm1, [r15 + pos_ *4]
	mulps xmm1, h_
	mulps xmm1, h_
	movdqu xmm3, [seis]
	divps xmm1, xmm3
	addps xmm0, xmm1
	mulps xmm0, wt_

	movdqu xmm1, [uno]
	subps xmm1, wt_

	mov r15, [mat_arr_ + offset_psix2_]
	movdqu xmm2, [r15 + pos_*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	mov r15, [mat_arr_ + offset_psix2_]
	mov r14, pos_
	movdqu [r15 + r14*4], xmm0


; --------------------Calculo de eje y----------------------------------
	
	call calcular_v 
	
	;Segundo calculo de delta: (1 + q * [V2.at(i, j - 1, k) - V2.at(i, j + 1, k)] + 6 * r)
	mov rax, pos_
	sub rax, offset_j_ ; aca tenemos pos para i, j-1, k
	mov r15, [mat_arr_ + offset_V2_]
	movdqu xmm0, [r15 + rax * 4]
	mov rax, pos_
	add rax, offset_j_ ; aca tenemos pos para i, j+1, k
	movdqu xmm1, [r15 + rax * 4]
	subps xmm0, xmm1
	mulps xmm0, q_
	movdqu xmm2, [uno]
	addps xmm0, xmm2
	movdqu xmm2, [seis]
	mulps xmm2, r_
	addps xmm0, xmm2 
	movdqu delta_, xmm0

	;from now on xmm0 acumulates the result
	xorps xmm0, xmm0

	;p1 = -U2.at(i + 1, j, k) + r
	mov r15, [mat_arr_ + offset_U2_]
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1
	subps xmm1, xmm2
	addps xmm1, r_
	;p1 * omy2.at(i + 1, j, k)
	mov r14, [mat_arr_ + offset_omy2_]
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p3 = U2.at(i - 1, j, k) + r
	mov rax, pos_
	sub rax, offset_i_				; rax = i-1,j,k
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p3 * omy2.at(i - 1, j, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p2 = -V2.at(i, j + 1, k) + r
	mov r15, [mat_arr_ + offset_V2_]	;r15 == V2
	mov rax, pos_
	add rax, offset_j_					; rax = i,j+1,k
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1
	subps xmm1, xmm2
	addps xmm1, r_
	;p2 * omy2.at(i, j + 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p4 = V2.at(i, j - 1, k) + r
	mov rax, pos_
	sub rax, offset_j_					; rax = i, j - 1, k
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p4 * omy2.at(i, j - 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p5 = W2.at(i, j, k + 1) + r
	mov r15, [mat_arr_ + offset_W2_]	;r15 == W2
	mov rax, pos_
	add rax, offset_k_					; rax = i, j, k + 1
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p5 * omy2.at(i, j, k + 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p6 = W2.at(i, j, k - 1) + r
	mov rax, pos_
	sub rax, offset_k_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p6 * omy2.at(i, j, k - 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1
 

	;(1.0 / q)*omy1.at(i, j, k)
	mov r15, [mat_arr_ + offset_omy1_]
	movdqu xmm1, [r15 + pos_*4]
	divps xmm1, q_
	addps xmm0, xmm1

	;ax1 = omx2.at(i, j, k) * V2.at(i + 1, j, k) 
	mov r15, [mat_arr_ + offset_omx2_]
	movdqu xmm1, [r15 + pos_*4]
	mov r14, [mat_arr_ + offset_V2_]
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

    ;ax2 = -omx2.at(i, j, k) * V2.at(i - 1, j, k)
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;ax3 = omz2.at(i, j, k) * V2.at(i, j, k + 1)
    mov r15, [mat_arr_ + offset_omz2_]
	movdqu xmm1, [r15 + pos_*4]
	mov rax, pos_
	add rax, offset_k_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

    ;ax4 = -omz2.at(i, j, k) * V2.at(i, j, k - 1)
	mov rax, pos_
	sub rax, offset_k_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;* wt * (q / delta)
	mulps xmm0, q_
	divps xmm0, delta_
	mulps xmm0, wt_

	;(1.0 - wt)*omy2(i,j,k)_old 
    mov r15, [mat_arr_ + offset_omy2_]
	movdqu xmm1, [r15 + pos_*4]
	movdqu xmm2, [uno]
	subps xmm2, wt_
	mulps xmm1, xmm2
	addps xmm0, xmm1

    mov r15, [mat_arr_ + offset_omy2_]
	mov r14, pos_
	shl r14, 2
	add r15, r14
	movdqu [r15], xmm0


;------------------------Eliptica eje y----------------------------
	mov r15, [mat_arr_ + offset_psiy2_]
	mov r14, pos_
	add r14, offset_i_
	movdqu xmm3,  [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	add r14, offset_j_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_i_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_j_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	add r14, offset_k_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_k_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r15, [mat_arr_ + offset_omy2_]
	movdqu xmm1, [r15 + pos_ *4]
	mulps xmm1, h_
	mulps xmm1, h_
	movdqu xmm3, [seis]
	divps xmm1, xmm3
	addps xmm0, xmm1
	mulps xmm0, wt_

	movdqu xmm1, [uno]
	subps xmm1, wt_
	mov r15, [mat_arr_ + offset_psiy2_]
	movdqu xmm2, [r15 + pos_*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	mov r15, [mat_arr_ + offset_psiy2_]
	mov r14, pos_
	movdqu [r15 + r14*4], xmm0


; --------------------Calculo de eje z----------------------------------
	
	call calcular_v 
	
	;Tercer calculo de delta: (1 + q * [W2.at(i, j, k - 1) - W2.at(i, j, k + 1)] + 6 * r)
	mov rax, pos_
	sub rax, offset_k_
	mov r15, [mat_arr_ + offset_W2_]
	movdqu xmm0, [r15 + rax * 4]
	mov rax, pos_
	add rax, offset_k_
	movdqu xmm1, [r15 + rax * 4]
	subps xmm0, xmm1
	mulps xmm0, q_
	movdqu xmm2, [uno]
	addps xmm0, xmm2
	movdqu xmm2, [seis]
	mulps xmm2, r_
	addps xmm0, xmm2 
	movdqu delta_, xmm0

 	;from now on xmm0 acumulates the result
 	xorps xmm0, xmm0

	;p1 = -U2.at(i + 1, j, k) + r
	mov r15, [mat_arr_ + offset_U2_]
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1
	subps xmm1, xmm2
	addps xmm1, r_
	;p1 * omz2.at(i + 1, j, k)
	mov r14, [mat_arr_ + offset_omz2_]
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p3 = U2.at(i - 1, j, k) + r
	mov rax, pos_
	sub rax, offset_i_				; rax = i-1,j,k
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p3 * omz2.at(i - 1, j, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p2 = -V2.at(i, j + 1, k) + r
	mov r15, [mat_arr_ + offset_V2_]	;r15 == V2
	mov rax, pos_
	add rax, offset_j_					; rax = i,j+1,k
	movdqu xmm2, [r15 + rax*4]
	xorps xmm1, xmm1
	subps xmm1, xmm2
	addps xmm1, r_
	;p2 * omz2.at(i, j + 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p4 = V2.at(i, j - 1, k) + r
	mov rax, pos_
	sub rax, offset_j_					; rax = i, j - 1, k
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p4 * omz2.at(i, j - 1, k)
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p5 = W2.at(i, j, k + 1) + r
	mov r15, [mat_arr_ + offset_W2_]	;r15 == W2
	mov rax, pos_
	add rax, offset_k_					; rax = i, j, k + 1
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p5 * omz2.at(i, j, k + 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	;p6 = W2.at(i, j, k - 1) + r
	mov rax, pos_
	sub rax, offset_k_
	movdqu xmm1, [r15 + rax*4]
	addps xmm1, r_
	;p6 * omz2.at(i, j, k - 1)
 	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1
 

	;(1.0 / q)*omz1.at(i, j, k)
	mov r15, [mat_arr_ + offset_omz1_]
	movdqu xmm1, [r15 + pos_*4]
	divps xmm1, q_
	addps xmm0, xmm1

	;ax1 = omx2.at(i, j, k) * W2.at(i + 1, j, k)
	mov r15, [mat_arr_ + offset_omx2_]
	movdqu xmm1, [r15 + pos_*4]
	mov r14, [mat_arr_ + offset_W2_]
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

	;ax2 = -omx2.at(i, j, k) * W2.at(i - 1, j, k)
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;ax3 = omy2.at(i, j, k) * W2.at(i, j + 1, k)
	mov r15, [mat_arr_ + offset_omy2_]
	movdqu xmm1, [r15 + pos_*4]
	mov rax, pos_
	add rax, offset_j_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	addps xmm0, xmm2

	;ax4 = -omy2.at(i, j, k) * W2.at(i, j - 1, k)
	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm2, [r14 + rax*4]
	mulps xmm2, xmm1
	subps xmm0, xmm2

	;* wt * (q / delta)
	mulps xmm0, q_
	divps xmm0, delta_
	mulps xmm0, wt_

	;(1.0 - wt)*omz2(i,j,k)_old 
	mov r15, [mat_arr_ + offset_omz2_]
	movdqu xmm1, [r15 + pos_*4]
	movdqu xmm2, [uno]
	subps xmm2, wt_
	mulps xmm1, xmm2
	addps xmm0, xmm1

	mov r15, [mat_arr_ + offset_omz2_]
	mov r14, pos_
	shl r14, 2
	add r15, r14
	movdqu [r15], xmm0


;------------------------Eliptica eje z----------------------------

	mov r15, [mat_arr_ + offset_psiz2_]
	mov r14, pos_
	add r14, offset_i_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	add r14, offset_j_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_i_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_j_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	add r14, offset_k_
	movdqu xmm3, [r15 + r14*4]
	addps xmm0, xmm3

	mov r14, pos_
	sub r14, offset_k_
	movdqu xmm3, [r15 + r14*4] 
	addps xmm0, xmm3

	mov r15, [mat_arr_ + offset_omz2_]
	movdqu xmm1, [r15 + pos_ *4]
	mulps xmm1, h_
	mulps xmm1, h_
	movdqu xmm3, [seis]
	divps xmm1, xmm3
	addps xmm0, xmm1
	mulps xmm0, wt_

	movdqu xmm1, [uno]
	subps xmm1, wt_
	mov r15, [mat_arr_ + offset_psiz2_]
	movdqu xmm2, [r15 + pos_*4]
	mulps xmm1, xmm2
	addps xmm0, xmm1

	mov r15, [mat_arr_ + offset_psiz2_]
	mov r14, pos_
	movdqu [r15 + r14*4], xmm0

	call calcular_v 

	mov rax, rdi

	add rsp, 8
	pop rbx
	pop r12
	pop r13
	pop r14
	pop r15
	pop rbp
	ret

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
calcular_v:
  	;int calcular_v (
    ;float* data,     rdi
    ;int    pos       rsi
    ;)

	push rbp
	mov rbp, rsp
	push r15  ;estos son los 5 registros a preservar... si es que los pisamos
	push r14
	push r13
	push r12
	push rbx
	sub rsp, 8  ;la pila debe estar alineada a 16 bits

	;---------------set U2------------------------
	mov r15, [mat_arr_ + offset_psiz2_]
	mov r14, pos_
	add r14, offset_j_
	movdqu xmm0, [r15+r14*4] ;(psiz2.at(i, j + 1, k)

	mov r14, pos_
	sub r14, offset_j_
	movdqu xmm3, [r15+r14*4]
	subps xmm0, xmm3 ; -psiz2.at(i, j - 1, k)

	mov r15, [mat_arr_ + offset_psiy2_]
	mov r14, pos_
	add r14, offset_k_
	movdqu xmm3, [r15+r14*4]
	subps xmm0, xmm3 ; -psiy2.at(i, j, k + 1)

	mov r14, pos_
	sub r14, offset_k_
	movdqu xmm3, [r15+r14*4] 
	addps xmm0, xmm3 ; -psiy2.at(i, j, k - 1)

	movdqu xmm3, [dos]
	divps xmm0, xmm3
	divps xmm0, h_

	mov r15, [mat_arr_ + offset_U2_]
	movdqu [r15 + pos_*4], xmm0

	;---------------set V2------------------------
	mov r15, [mat_arr_ + offset_psix2_]
	mov r14, pos_
	add r14, offset_k_
	movdqu xmm0, [r15+r14*4] ;psix2.at(i, j, k + 1)	

	mov r14, pos_
	sub r14, offset_k_
	movdqu xmm3, [r15+r14*4] 
	subps xmm0, xmm3 ;psix2.at(i, j, k - 1)	

	mov r15, [mat_arr_ + offset_psiz2_]
	mov r14, pos_
	add r14, offset_i_
	movdqu xmm3, [r15+r14*4]
	subps xmm0, xmm3 ;psiz2.at(i+1, j, k)	

	mov r14, pos_
	sub r14, offset_i_
	movdqu xmm3, [r15+r14*4]
	addps xmm0, xmm3 ;psiz2.at(i-1, j, k)

	movdqu xmm3, [dos]
	divps xmm0, xmm3
	divps xmm0, h_

	mov r15, [mat_arr_ + offset_V2_]
	movdqu [r15 + pos_*4], xmm0

	;---------------set W2------------------------
	mov r15, [mat_arr_ + offset_psiy2_]
	mov r14, pos_
	add r14, offset_i_
	movdqu xmm0, [r15+r14*4] ;psiy2.at(i + 1, j, k)

	mov r14, pos_
	sub r14, offset_i_
	movdqu xmm3, [r15+r14*4]
	subps xmm0, xmm3 ;psiy2.at(i - 1, j, k)

	mov r15, [mat_arr_ + offset_psix2_]
	mov r14, pos_
	add r14, offset_j_
	movdqu xmm3, [r15+r14*4]
	subps xmm0, xmm3 ;psix2.at(i, j + 1, k)

	mov r14, pos_
	sub r14, offset_j_
	movdqu xmm3, [r15+r14*4]
	addps xmm0, xmm3 ;psix2.at(i, j + 1, k)

	movdqu xmm3, [dos]
	divps xmm0, xmm3
	divps xmm0, h_

	mov r15, [mat_arr_ + offset_W2_]
	movdqu [r15 + pos_*4], xmm0

    add rsp, 8
	pop rbx
	pop r12
	pop r13
	pop r14
	pop r15
	pop rbp
	ret
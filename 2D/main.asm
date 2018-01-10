[BITS 64]

extern printf
global vvp_asm


section .data
copy_float_value_mask: DB 0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3
eltest: DB 0x05
uno: DD 1.0, 1.0, 1.0, 1.0
dos: DD 2.0, 2.0, 2.0, 2.0

section .text
%define dt_ xmm0
%define dx_ xmm1
%define dy_ xmm2
%define rho_ xmm3
%define nu_ xmm4
%define u1ij_ xmm5
%define v1ij_ xmm6

%define u1p_ r15
%define v1p_ r14
%define offset_i_ rdx
%define offset_j_ 1
%define mat_arr_ rdi
%define pos_ rsi

%define offset_U2_ 0
%define offset_V2_ 8
%define offset_P2_ 16
%define offset_U1_ 24
%define offset_V1_ 32
%define offset_P1_ 40


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

vvp_asm:	
	push rbp
	mov rbp, rsp
	push r15  ;estos son los 5 registros a preservar... si es que los pisamos
	push r14
	push rax
	push rbx
	;sub rsp, 8  ;la pila debe estar alineada a 16 bits

  ;int vvp_asm  (
    ;float* data,     rdi
    ;int    pos       rsi
    ;int 	offset i  rdx
    ;float  dt		  xmm0
    ;float  dx		  xmm1
    ;float  dy		  xmm2
    ;float  rho		  xmm3
    ;float  nu		  xmm4
  ;)


	;Extending constants to the whole register. 
	movdqu xmm15, [copy_float_value_mask]
	pshufb dt_, xmm15
	pshufb dx_, xmm15
	pshufb dy_, xmm15

	;setting pointers
	mov u1p_, [mat_arr_ + offset_U1_]
	movdqu u1ij_, [u1p_ + pos_]

	mov v1p_, [mat_arr_ + offset_V1_]
	movdqu v1ij_, [v1p_ + pos_]

	;-------------Setting U2-------------------;
	;from now on xmm15 acumulates the result
	xorps xmm15, xmm15

	;U1.at(i,j)
	movdqu xmm15, u1ij_ ;176

	;- U1.at(i, j) *(dt/dx)
	xorps xmm14, xmm14
	subps xmm14, u1ij_
	mulps  xmm14, dt_
	divps xmm14, dx_

	;U1.at(i, j)-U1.at(i - 1, j)
	movdqu xmm13, u1ij_
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm12, [u1p_ + rax]
	subps xmm13, xmm12

	mulps xmm14, xmm13
	addps xmm15, xmm14	;hasta 177

	;- V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1)) 
	movdqu xmm14, v1ij_
	mulps xmm14, dt_
	divps xmm14, dy_

	movdqu xmm13, u1ij_
	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm12, [u1p_ + rax]
	subps xmm13, xmm12
	mulps xmm14, xmm13

	subps xmm15, xmm14  ;hasta 178

	;(dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
	movdqu xmm14, dt_
	movdqu xmm13, rho_
	mulps xmm13, [dos]
	mulps xmm13, dx_
	divps xmm14, xmm13

	mov r13, [mat_arr_ + offset_P1_] 
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm13, [r13 + rax]
	sub rax, offset_i_
	sub rax, offset_i_
	movdqu xmm12, [r13 + rax]
	subps xmm13, xmm12

	mulps xmm14, xmm13	;hasta 179 
	subps xmm15, xmm14
;+ nu * (
;		(dt / (dx * dx)) * (U1.at(i + 1, j)
;- 2 * U1.at(i, j) + U1.at(i - 1, j)) 

	
	mov rax, pos_
	add rax, offset_i_
	movdqu xmm14, [u1p_ + rax]
	subps xmm14, u1ij_
	subps xmm14, u1ij_

	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm13, [u1p_ + rax]
	addps xmm14, xmm13

	mulps xmm14, dt_
	divps xmm14, dx_
	divps xmm14, dx_

;  + (dt / (dy * dy)) * (U1.at(i, j + 1)
; 	- 2 * U1.at(i, j) + U1.at(i, j - 1)) )

	mov rax, pos_
	add rax, offset_j_
	movdqu xmm12, [u1p_ + rax]
	subps xmm12, u1ij_
	subps xmm12, u1ij_

	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm13, [u1p_ + rax]
	addps xmm12, xmm13

	mulps xmm12, dt_
	divps xmm12, dy_
	divps xmm12, dy_

	;adding the two terms
	addps xmm14, xmm12
	mulps xmm14, nu_ 
	addps xmm15, xmm14 ;181

	mov r13, [mat_arr_ + offset_U2_]
	add r13, pos_
	movdqu [r13], xmm15 ;182

	;-------------Setting V2-------------------;
	;from now on xmm15 acumulates the result
	xorps xmm15, xmm15

	;U1.at(i,j)
	movdqu xmm15, v1ij_ ;184

	;- U1.at(i, j) *(dt/dx)
	xorps xmm14, xmm14
	subps xmm14, u1ij_
	mulps  xmm14, dt_
	divps xmm14, dx_

	;V1.at(i, j)-V1.at(i - 1, j)
	movdqu xmm13, v1ij_
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm12, [v1p_ + rax]
	subps xmm13, xmm12

	mulps xmm14, xmm13
	addps xmm15, xmm14	;hasta 185

	;- V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i-1, j)) 
	movdqu xmm14, v1ij_
	mulps xmm14, dt_
	divps xmm14, dy_

	movdqu xmm13, v1ij_
	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm12, [v1p_ + rax]
	subps xmm13, xmm12
	mulps xmm14, xmm13

	subps xmm15, xmm14  ;hasta 186

	;(dt / (rho * 2 * dy)) * (P1.at(i, j+1) - P1.at(i, j-1))
	movdqu xmm14, dt_
	movdqu xmm13, rho_
	mulps xmm13, [dos]
	mulps xmm13, dy_
	divps xmm14, xmm13

	mov r13, [mat_arr_ + offset_P1_] 
	mov rax, pos_
	add rax, offset_j_
	movdqu xmm13, [r13 + rax]
	sub rax, offset_j_
	sub rax, offset_j_
	movdqu xmm12, [r13 + rax]
	subps xmm13, xmm12

	mulps xmm14, xmm13	;hasta 187 
	subps xmm15, xmm14

;+ nu * (
;		(dt / (dx * dx)) * (V1.at(i + 1, j)
;- 2 * V1.at(i, j) + U1.at(i - 1, j)) 

	mov rax, pos_
	add rax, offset_i_
	movdqu xmm14, [v1p_ + rax]
	subps xmm14, v1ij_
	subps xmm14, v1ij_

	mov rax, pos_
	sub rax, offset_i_
	movdqu xmm13, [v1p_ + rax]
	addps xmm14, xmm13

	mulps xmm14, dt_
	divps xmm14, dx_
	divps xmm14, dx_

;  + (dt / (dy * dy)) * (V1.at(i, j + 1)
; 	- 2 * V1.at(i, j) + V1.at(i, j - 1)) )

	mov rax, pos_
	add rax, offset_j_
	movdqu xmm12, [v1p_ + rax]
	subps xmm12, v1ij_
	subps xmm12, v1ij_

	mov rax, pos_
	sub rax, offset_j_
	movdqu xmm13, [v1p_ + rax]
	addps xmm12, xmm13

	mulps xmm12, dt_
	divps xmm12, dy_
	divps xmm12, dy_

	;adding the two terms
	addps xmm14, xmm12
	mulps xmm14, nu_ 
	addps xmm15, xmm14 ;189

	mov r13, [mat_arr_ + offset_V2_]
	add r13, pos_
	movdqu [r13], xmm15 ;190

	;-------------Setting Final Term-------------------;











	;add rsp, 8
	pop rbx
	pop rax
	pop r14
	pop r15
	pop rbp
	ret

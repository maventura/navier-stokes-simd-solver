[BITS 64]

extern printf
global vvp_asm


section .data
copy_float_value_mask: DB 0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3
eltest: DB 0x05
wt: DD 0.8
uno: DD 1.0, 1.0, 1.0, 1.0
seis: DD 6.0, 6.0, 6.0, 6.0

section .text
%define dos 0x4000000000000000
%define r_ xmm10
%define q_ xmm11
%define h_ xmm12
%define wt_ xmm13
%define delta_ xmm14
%define qinv_ xmm15
%define offset_i_ rdx
%define offset_j_ rcx
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
	movdqu h_, xmm0
	movdqu xmm0, [copy_float_value_mask]
	pshufb h_, xmm0

	xorps q_, q_
	movdqu q_, xmm0
	movdqu xmm0, [copy_float_value_mask]
	pshufb q_, xmm0

	xorps wt_, wt_
	movdqu wt_, [wt]
	movdqu xmm0, [copy_float_value_mask]
	pshufb wt_, xmm0


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

	movdqu qinv_, [uno]
	divps qinv_, q_


	;from now on xmm0 acumulates the result
	xorps xmm0, xmm0
	
	;P1)
	;xmm0 //res/acum p_i/omx2(i, j, k)    	double p1 = (-U2.at(i + 1, j, k) + r);
	mov r15, [mat_arr_ + offset_U2_]	;r15 == U2
	mov rax, pos_
	add rax, offset_i_					; rax = i+1,j,k
	movdqu xmm2, [r15 + rax*4]
	subps xmm1, xmm2

	addps xmm1, r_

	mov r14, [mat_arr_ + offset_omx2_]	;r14 = omx2
	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	movdqu xmm0, xmm1

	;P3)
	;double p3 = (U2.at(i - 1, j, k) + r);   

	;r15 == U2
	mov rax, pos_
	sub rax, offset_i_				; rax = i-1,j,k
	movdqu xmm1, [r15 + rax*4]

	addps xmm1, r_

	movdqu xmm2, [r14 + rax*4]
	mulps xmm1, xmm2
	movdqu xmm0, xmm1





	; sub xmm1, p_i ;(las primeras dos veces)    	double p2 = (-V2.at(i, j + 1, k) + r)	
	; add xmm1, p_i ;(las 4 veces restantes)	    										
	; add xmm1, r_(r)	;                                            
	; mov xmm2, om_i    ;	double p3 = (U2.at(i - 1, j, k) + r);   
	; mul xmm1, xmm2    ;	double p4 = (V2.at(i, j - 1, k) + r);   
	; add xmm0, xmm1   ; 	double p5 = (W2.at(i, j, k + 1) + r);   
	; ;"c" = xmm internal index.    	double p6 = (W2.at(i, j, k - 1) + r);   
 ;                                            ;no tocar xmm0

	; xor xmm1, xmm1 ; axs/acum ax_i/res
    
 ;    mov xmm2, omy3
	; mov xmm4, U2
	; mul xmm2, xmm4
	; add xmm1, xmm2 ;double ax1 = omy2.at(i, j, k) * U2.at(i, j + 1, k); 
 ;    ;double ax2 = -omy2.at(i, j, k) * U2.at(i, j - 1, k);//same
 ;    ;double ax3 = omz2.at(i, j, k) * U2.at(i, j, k + 1);//same
 ;    ;double ax4 = -omz2.at(i, j, k) * U2.at(i, j, k - 1);//same

 ;    ;double axs = ax1 + ax2 + ax3 + ax4;

 ;    ;omx2.set(i, j, k, p1 * omx2.at(i + 1, j, k) + p2 * omx2.at(i, j + 1, k)//xmm0
 ;             + p3 * omx2.at(i - 1, j, k) + p4 * omx2.at(i, j - 1, k)//xmm0
 ;             + p5 * omx2.at(i, j, k + 1) + p6 * omx2.at(i, j, k - 1)//xmm0
 ;             + (1.0 / q)*omx1.at(i, j, k)   
 ;    mov qinv_, 1/q ;(FIJO)
 ;    mov xmm2, omx1
 ;    mul xmm2, qinv_
 ;             ;+ axs);//xmm1

    
 ;    omx2.set(i, j, k, omx2.at(i, j, k) * (q / delta)); //add xmm0, xmm1
 ;                                                        //add xmm0, xmm2    //  xmm0 = set. (xmm1/3 libres)
 ;                                                        //div xmm0, qinv_   //  /= delta
 ;                                                        //div xmm0, delta_   //  *=(1/q)

 ;        omx2.set(i, j, k, (1.0 - wt)*aux_omx2 + wt * omx2.at(i, j, k));//mov wt_, wt (=0.8 FIJO)
 ;                                                                      //mul xmm0, wt_
 ;                                                                      //mov xmm1, 0x0001
 ;                                                                      //sub xmm1, wt_
 ;                                                                      //mov xmm2, omx2(i,j,k)
 ;                                                                      //mul xmm2, xmm1
 ;                                                                      //add xmm0, xmm2 

 ;                                                            //mov [omx2], xmm0


 ;    //Eliptica eje x.
 ;    //p1=p2=p3=p4=p5=p6=1/6.0;                                    
 ;    psix2.set(i, j, k, (psix2.at(i + 1, j, k) + psix2.at(i, j + 1, k)   //mov xmm0, psix2 (levanto uno, una vez)
 ;                        + psix2.at(i - 1, j, k) + psix2.at(i, j - 1, k) //mov xmm1, pxix2 (levanto otro, aca itera)
 ;                        + psix2.at(i, j, k + 1) + psix2.at(i, j, k - 1) //add xmm0, xmm1 (sumo y vuelvo arriba)
 ;                        + h * h * omx2.at(i, j, k)) / 6.0);             // mov h_, h (DEPENDE PARAMS DE ENTRADA, pero nunca cambia)

	; mov xmm1, h_
	; mul xmm1, xmmm12
	; mov xmm2, omx2
	; mul xmm1, xmm2
	; mov xmmverde, 6.0 4 veces
	; div xmm1, xmmverde 
	; add xmm0, xmm1

 ;    psix2.set(i, j, k, (1.0 - wt)*aux_psix2 + wt * psix2.at(i, j, k));  //mul xmm0, wt_
 ;                                                                      //mov xmm1, 0x0001
 ;                                                                      //sub xmm1, wt_
 ;                                                                      //mov xmm2, psix2(i,j,k)
 ;                                                                      //mul xmm2, xmm1
 ;                                                                      //add xmm0, xmm2 
 ;                                                                     //mov [psix2], xmm0


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

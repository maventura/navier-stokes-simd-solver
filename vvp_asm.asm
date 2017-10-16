[BITS 64]

extern printf
global vvp_asm


section .data
copy_float_value_mask: DB 0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3,0x0, 0x1, 0x2, 0x3
eltest: DB 0x05
wt: DD 0.8
uno: DD 1.0, 1.0, 1.0, 1.0

section .text
%define dos 0x4000000000000000

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

	;en xmm0 esta el offset

	xorps xmm10, xmm10
	movdqu xmm10, xmm0
	movdqu xmm0, [copy_float_value_mask]  ; r esta 4 veces en xmm10
	pshufb xmm10, xmm0

	xorps xmm12, xmm12
	movdqu xmm12, xmm1
	movdqu xmm1, [copy_float_value_mask]  ; h o dx esta 4 veces en xmm11
	pshufb xmm12, xmm1

	xorps xmm11, xmm11
	movdqu xmm11, xmm2
	movdqu xmm2, [copy_float_value_mask]  ; q esta 4 veces en xmm11
	pshufb xmm11, xmm2

	xorps xmm13, xmm13
	movdqu xmm13, [wt]
	movdqu xmm1, [copy_float_value_mask]  ; wt o 0.8 esta 4 veces en xmm13
	pshufb xmm13, xmm1


	;calculo de delta: (1 + q * [U2.at(i - 1, j, k) - U2.at(i + 1, j, k)] + 6 * r)
	mov rax, rsi
	sub rax, rcx ; aca tenemos pos para i-1, j, k
	movdqu xmm1, [rdi + rax * 4]
	mov rax, rsi
	add rax, rcx ; aca tenemos pos para i+1, j, k
	movdqu xmm2, [rdi + rax * 4]

	subps xmm1, xmm2
	mulps xmm1, xmm11
	addps xmm1, [uno]


	; mov xmm14, delta
	; xor xmm1, xmm1
	; ;xmm1 //res/acum p_i/omx2(i, j, k)    	double p1 = (-U2.at(i + 1, j, k) + r);  
	; xor xmm2, xmm2    																					
	; sub xmm2, p_i ;(las primeras dos veces)    	double p2 = (-V2.at(i, j + 1, k) + r)	
	; add xmm2, p_i ;(las 4 veces restantes)	    										
	; add xmm2, xmm10(r)	;                                            
	; mov xmm3, om_i    ;	double p3 = (U2.at(i - 1, j, k) + r);   
	; mul xmm2, xmm3    ;	double p4 = (V2.at(i, j - 1, k) + r);   
	; add xmm1, xmm2   ; 	double p5 = (W2.at(i, j, k + 1) + r);   
	; ;"c" = xmm internal index.    	double p6 = (W2.at(i, j, k - 1) + r);   
 ;                                            ;no tocar xmm1

	; xor xmm2, xmm2 ; axs/acum ax_i/res
    
 ;    mov xmm3, omy3
	; mov xmm4, U2
	; mul xmm3, xmm4
	; add xmm2, xmm3 ;double ax1 = omy2.at(i, j, k) * U2.at(i, j + 1, k); 
 ;    ;double ax2 = -omy2.at(i, j, k) * U2.at(i, j - 1, k);//same
 ;    ;double ax3 = omz2.at(i, j, k) * U2.at(i, j, k + 1);//same
 ;    ;double ax4 = -omz2.at(i, j, k) * U2.at(i, j, k - 1);//same

 ;    ;double axs = ax1 + ax2 + ax3 + ax4;

 ;    ;omx2.set(i, j, k, p1 * omx2.at(i + 1, j, k) + p2 * omx2.at(i, j + 1, k)//xmm1
 ;             + p3 * omx2.at(i - 1, j, k) + p4 * omx2.at(i, j - 1, k)//xmm1
 ;             + p5 * omx2.at(i, j, k + 1) + p6 * omx2.at(i, j, k - 1)//xmm1
 ;             + (1.0 / q)*omx1.at(i, j, k)   
 ;    mov xmm15, 1/q ;(FIJO)
 ;    mov xmm3, omx1
 ;    mul xmm3, xmm15
 ;             ;+ axs);//xmm2

    
 ;    omx2.set(i, j, k, omx2.at(i, j, k) * (q / delta)); //add xmm1, xmm2
 ;                                                        //add xmm1, xmm3    //  xmm1 = set. (xmm2/3 libres)
 ;                                                        //div xmm1, xmm15   //  /= delta
 ;                                                        //div xmm1, xmm14   //  *=(1/q)

 ;        omx2.set(i, j, k, (1.0 - wt)*aux_omx2 + wt * omx2.at(i, j, k));//mov xmm13, wt (=0.8 FIJO)
 ;                                                                      //mul xmm1, xmm13
 ;                                                                      //mov xmm2, 0x0001
 ;                                                                      //sub xmm2, xmm13
 ;                                                                      //mov xmm3, omx2(i,j,k)
 ;                                                                      //mul xmm3, xmm2
 ;                                                                      //add xmm1, xmm3 

 ;                                                            //mov [omx2], xmm1


 ;    //Eliptica eje x.
 ;    //p1=p2=p3=p4=p5=p6=1/6.0;                                    
 ;    psix2.set(i, j, k, (psix2.at(i + 1, j, k) + psix2.at(i, j + 1, k)   //mov xmm1, psix2 (levanto uno, una vez)
 ;                        + psix2.at(i - 1, j, k) + psix2.at(i, j - 1, k) //mov xmm2, pxix2 (levanto otro, aca itera)
 ;                        + psix2.at(i, j, k + 1) + psix2.at(i, j, k - 1) //add xmm1, xmm2 (sumo y vuelvo arriba)
 ;                        + h * h * omx2.at(i, j, k)) / 6.0);             // mov xmm12, h (DEPENDE PARAMS DE ENTRADA, pero nunca cambia)

	; mov xmm2, xmm12
	; mul xmm2, xmmm12
	; mov xmm3, omx2
	; mul xmm2, xmm3
	; mov xmmverde, 6.0 4 veces
	; div xmm2, xmmverde 
	; add xmm1, xmm2

 ;    psix2.set(i, j, k, (1.0 - wt)*aux_psix2 + wt * psix2.at(i, j, k));  //mul xmm1, xmm13
 ;                                                                      //mov xmm2, 0x0001
 ;                                                                      //sub xmm2, xmm13
 ;                                                                      //mov xmm3, psix2(i,j,k)
 ;                                                                      //mul xmm3, xmm2
 ;                                                                      //add xmm1, xmm3 

 ;                                                                     //mov [psix2], xmm1






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;testeo de memoria (no borrar, a martin le duele);;;;;;;;;;;;;;;;;;;;;;;;;;;;
;rdi pointer to data.
;rsi x
;rdx y
;rcx z

;mov r9, 0x4000000000000000
;mov [rdi], r9

;mov dword [rdi], 0
;mov dword [rdi+4], 0x403F0000

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

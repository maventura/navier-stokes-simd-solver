[BITS 64]

global vvp_asm

%define dos 0x4000000000000000

section .text

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vvp_asm:	

	push rbp
	mov rbp, rsp
	push r8
	push r9
	push r10

	;chequear si hace falta asegurar los xmmi
	;en xmm0 esta el offset
	;en  rdi el puntero a la matriz

	movd xmm15, 0.8
	mov xmm10, ; r = dt / (Re * h * h)


	mov xmm14, delta
	xor xmm1, xmm1
    ;xmm10,  r (FIJO) pasar a xmm15
    double p1 = (-U2.at(i + 1, j, k) + r);  //xmm1 //res/acum p_i/omx2(i, j, k)
											;xor xmm2, xmm2    									
    double p2 = (-V2.at(i, j + 1, k) + r)	;sub xmm2, p_i (las primeras dos veces)
    										;add xmm2, p_i (las 4 veces restantes)
                                            ;add xmm2, xmm10(r)
    double p3 = (U2.at(i - 1, j, k) + r);   //mov xmm3, om_i
    double p4 = (V2.at(i, j - 1, k) + r);   //mul xmm2, xmm3
    double p5 = (W2.at(i, j, k + 1) + r);   //add xmm1, xmm2
    double p6 = (W2.at(i, j, k - 1) + r);   //"c" = xmm internal index.
                                            //no tocar xmm1

                                                        //xor xmm2, xmm2 // axs/acum ax_i/res
    double ax1 = omy2.at(i, j, k) * U2.at(i, j + 1, k); //mov xmm3, omy3
                                                        //mov xmm4, U2
                                                        //mul xmm3, xmm4
                                                        //add xmm2, xmm3
    double ax2 = -omy2.at(i, j, k) * U2.at(i, j - 1, k);//same
    double ax3 = omz2.at(i, j, k) * U2.at(i, j, k + 1);//same
    double ax4 = -omz2.at(i, j, k) * U2.at(i, j, k - 1);//same

    double axs = ax1 + ax2 + ax3 + ax4;

    omx2.set(i, j, k, p1 * omx2.at(i + 1, j, k) + p2 * omx2.at(i, j + 1, k)//xmm1
             + p3 * omx2.at(i - 1, j, k) + p4 * omx2.at(i, j - 1, k)//xmm1
             + p5 * omx2.at(i, j, k + 1) + p6 * omx2.at(i, j, k - 1)//xmm1
             + (1.0 / q)*omx1.at(i, j, k)   //mov xmm15, 1/q (FIJO)
                                            //mov xmm3, omx1
                                            //mul xmm3, xmm15
             + axs);//xmm2

    
    omx2.set(i, j, k, omx2.at(i, j, k) * (q / delta)); //add xmm1, xmm2
                                                        //add xmm1, xmm3    //  xmm1 = set. (xmm2/3 libres)
                                                        //div xmm1, xmm15   //  /= delta
                                                        //div xmm1, xmm14   //  *=(1/q)

        omx2.set(i, j, k, (1.0 - wt)*aux_omx2 + wt * omx2.at(i, j, k));//mov xmm13, wt (=0.8 FIJO)
                                                                      //mul xmm1, xmm13
                                                                      //mov xmm2, 0x0001
                                                                      //sub xmm2, xmm13
                                                                      //mov xmm3, omx2(i,j,k)
                                                                      //mul xmm3, xmm2
                                                                      //add xmm1, xmm3 

                                                            //mov [omx2], xmm1


    //Eliptica eje x.
    //p1=p2=p3=p4=p5=p6=1/6.0;                                    
    psix2.set(i, j, k, (psix2.at(i + 1, j, k) + psix2.at(i, j + 1, k)   //mov xmm1, psix2 (levanto uno, una vez)
                        + psix2.at(i - 1, j, k) + psix2.at(i, j - 1, k) //mov xmm2, pxix2 (levanto otro, aca itera)
                        + psix2.at(i, j, k + 1) + psix2.at(i, j, k - 1) //add xmm1, xmm2 (sumo y vuelvo arriba)
                        + h * h * omx2.at(i, j, k)) / 6.0);             // mov xmm12, h (DEPENDE PARAMS DE ENTRADA, pero nunca cambia)

                                                                        //mov xmm2, xmm12
                                                                        //mul xmm2, xmmm12
                                                                        //mov xmm3, omx2
                                                                        //mul xmm2, xmm3
                                                                        //div xmm2, xmm11 (=6.0)
                                                                        //add xmm1, xmm2

    psix2.set(i, j, k, (1.0 - wt)*aux_psix2 + wt * psix2.at(i, j, k));  //mul xmm1, xmm13
                                                                      //mov xmm2, 0x0001
                                                                      //sub xmm2, xmm13
                                                                      //mov xmm3, psix2(i,j,k)
                                                                      //mul xmm3, xmm2
                                                                      //add xmm1, xmm3 

                                                                     //mov [psix2], xmm1






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
pop r10
pop r9
pop r8
pop rbp
ret
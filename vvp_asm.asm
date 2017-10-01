[BITS 64]

global vvp_asm

%define dos 0x4000000000000000

section .text
vvp_asm:

push rbp
mov rbp, rsp
push r8
push r10
push r9

;rdi pointer to data.
;rsi x
;rdx y
;rcx z

;mov r9, 0x4000000000000000
;mov [rdi], r9

mov dword [rdi], 0
mov dword [rdi+4], 0x403F0000

mov rax, rdi
pop r9
pop r10
pop r8
pop rbp
ret